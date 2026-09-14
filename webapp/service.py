"""Library bootstrap, NIM health checks, and a serial durable prediction worker."""
import csv
import hashlib
import json
import logging
import threading
import urllib.error
import urllib.request
from pathlib import Path
from webapp.domain import molecule, payload, scores
from webapp.store import now, uid
from webapp.benchmark import Benchmarks
from scripts.download_bindingdb import prepare_sample

log = logging.getLogger(__name__)


class Service:
    def __init__(self, store, endpoint, bootstrap=True):
        self.store = store
        self.endpoint = endpoint.rstrip('/')
        self.bootstrap_enabled = bootstrap
        self.stop = threading.Event()
        self.bootstrap_lock = threading.Lock()
        self.ready = False
        self.health_message = 'Connecting to Boltz-2'
        self.threads = []
        self.benchmarks = Benchmarks(store,self)

    def start(self, worker=True, prepare=True):
        self.is_worker=worker
        if worker:
            self.store.recover()
        if prepare:
            self.benchmarks.recover()
        for target in ([self.health_loop, self.work] if worker else [self.health_loop]):
            t=threading.Thread(target=target,daemon=True)
            self.threads.append(t)
            t.start()
        if self.bootstrap_enabled:
            self.bootstrap_async()

    def close(self):
        self.stop.set()
        for t in self.threads:
            t.join(timeout=1)

    def health_loop(self):
        while not self.stop.is_set():
            try:
                with urllib.request.urlopen(self.endpoint+'/v1/health/ready', timeout=4) as r:
                    self.ready = r.status == 200
                self.health_message = 'Boltz-2 ready' if self.ready else 'Model is starting'
            except OSError:
                self.ready = False
                self.health_message = 'Model is starting or unavailable'
            if getattr(self,'is_worker',False):
                self.store.set_meta('worker',dict(heartbeat=now(),ready=self.ready))
            self.stop.wait(5)

    def bootstrap_async(self):
        if not self.bootstrap_lock.acquire(blocking=False):
            return
        threading.Thread(target=self.bootstrap,daemon=True).start()

    def bootstrap(self):
        try:
            if self.store.meta('bootstrap',{}).get('status') == 'complete':
                return
            self.store.set_meta('bootstrap',dict(status='loading',message='Preparing five BindingDB example ligands'))
            folder = self.store.folder/'data'
            prepare_sample(folder)
            # Verify this release before importing records into the library.
            archive = folder/'raw'/'BindingDB_BindingDB_Articles_202609_tsv.zip'
            if hashlib.md5(archive.read_bytes()).hexdigest() != '8c41a1fcf8b828d99070f4bca6bd7a86':
                archive.unlink(missing_ok=True)
                raise ValueError('BindingDB checksum did not match the pinned release. Retry the download.')
            with (folder/'BindingDB_top5.tsv').open() as f:
                rows = list(csv.DictReader(f,delimiter='\t'))
            prepared = [(row,molecule(row['Ligand SMILES'])) for row in rows]
            with self.store.connect() as db:
                db.execute('INSERT OR IGNORE INTO targets VALUES (?,?,?,?)',
                           ('src','Human Src',rows[0]['BindingDB Target Chain Sequence 1'],'P12931'))
                for row,(smiles,svg,mw) in prepared:
                    id='bdb-'+row['BindingDB MonomerID']
                    db.execute('INSERT OR IGNORE INTO ligands VALUES (?,?,?,?,?)',
                               (id,row['BindingDB MonomerID'],smiles,svg,mw))
                    actual=db.execute('SELECT id FROM ligands WHERE smiles=?',(smiles,)).fetchone()['id']
                    db.execute('INSERT OR IGNORE INTO measurements VALUES (?,?,?,?,?,?)',
                               ('src',actual,'IC50',row['IC50 (nM)'].strip(),'nM',row['Article DOI']))
            self.store.set_meta('bootstrap',dict(status='complete',message='BindingDB example ready'))
        except Exception as e:
            log.exception('Example import failed')
            self.store.set_meta('bootstrap',dict(status='failed',message=str(e)[:400]))
        finally:
            self.bootstrap_lock.release()

    def work(self):
        while not self.stop.is_set():
            try:
                if self.ready and not self.store.meta('paused',False):
                    job = self.store.claim()
                    if job:
                        self.predict(job)
                        continue
            except Exception:
                log.exception('Queue worker error')
                self.store.set_meta('paused',True)
            self.stop.wait(1)

    def predict(self, job):
        folder=self.store.folder/'artifacts'/job['id']/str(job['attempt'])
        folder.mkdir(parents=True,exist_ok=True)
        snapshot=json.loads(job['snapshot'])
        request=payload(snapshot)
        (folder/'request.json').write_text(json.dumps(request,indent=2))
        try:
            req=urllib.request.Request(self.endpoint+'/biology/mit/boltz2/predict',
                data=json.dumps(request).encode(),headers={'Content-Type':'application/json'})
            with urllib.request.urlopen(req,timeout=7200) as r:
                raw=r.read()
            (folder/'response.json').write_bytes(raw)
            response=json.loads(raw)
            result=scores(response,snapshot['model_version'])
            structures=response.get('structures') or []
            if structures and isinstance(structures[0].get('structure'),str):
                (folder/'complex.cif').write_text(structures[0]['structure'])
            (folder/'scores.json').write_text(json.dumps(result,indent=2))
            self.store.finish(job,'succeeded',score=result)
        except urllib.error.HTTPError as e:
            detail=e.read().decode(errors='replace')[:8000]
            (folder/'error.txt').write_text(detail)
            self.store.finish(job,'failed',error=f'NIM returned HTTP {e.code}. {detail[:500]}',pause=e.code>=500)
        except (urllib.error.URLError,TimeoutError,ConnectionError,OSError) as e:
            self.store.finish(job,'interrupted',error='Connection ended before completion was confirmed. Restart the NIM before resuming the queue.',pause=True)
            (folder/'error.txt').write_text(str(e))
        except Exception as e:
            self.store.finish(job,'failed',error=str(e)[:600])
            (folder/'error.txt').write_text(str(e))
