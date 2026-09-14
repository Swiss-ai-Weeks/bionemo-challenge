"""Reproducible, streaming BindingDB sampling and persistent benchmark analysis."""
from collections import Counter, defaultdict
import csv
from datetime import datetime
from functools import lru_cache
import hashlib
import io
import json
import math
from pathlib import Path
import random
import re
import statistics
import threading
import zipfile
from rdkit import Chem
from webapp.domain import AA, DEFAULT_SETTINGS, MODEL_VERSION, CONVERSION
from webapp.store import now, uid
from scripts.download_bindingdb import ARCHIVE, MEMBER, prepare_sample

EXACT = re.compile(r'(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?')


@lru_cache(maxsize=50000)
def canonical_smiles(text):
    if not text or len(text)>2000:
        return None
    m=Chem.MolFromSmiles(text)
    if m is None or not 1<=m.GetNumHeavyAtoms()<=128 or len(Chem.GetMolFrags(m))!=1:
        return None
    return Chem.MolToSmiles(m)


def sample_rows(rows, count, seed, endpoints, min_residues, max_residues, progress=None, stop=None):
    """Uniform reservoir over distinct (sequence, canonical SMILES, endpoint).

    For repeated measurements of a selected pair, independently choose one exact
    label uniformly among that pair's eligible rows. Store selected records only.
    """
    rng=random.Random(seed)
    label_rng=random.Random(seed ^ 0xBDB)
    counts=Counter()
    seen={}
    selected={}
    slots=[]
    for row_number,row in enumerate(rows,1):
        counts['rows_scanned']+=1
        if row_number%1000==0:
            if stop and stop(): raise InterruptedError('Sampling cancelled.')
            if progress: progress(dict(counts))
        seq=(row.get('BindingDB Target Chain Sequence 1') or '').strip()
        if (row.get('Number of Protein Chains in Target (>1 implies a multichain complex)') or '').strip()!='1':
            counts['excluded_multichain']+=1;continue
        if not min_residues<=len(seq)<=max_residues or set(seq)-AA:
            counts['excluded_sequence']+=1;continue
        labels=[]
        for endpoint in endpoints:
            text=(row.get(endpoint+' (nM)') or '').strip()
            if not text: continue
            if not EXACT.fullmatch(text):
                counts['excluded_censored_or_invalid_labels']+=1;continue
            number=float(text)
            if not math.isfinite(number) or number<=0:
                counts['excluded_censored_or_invalid_labels']+=1;continue
            labels.append((endpoint,number,text))
        if not labels:
            counts['excluded_no_exact_label']+=1;continue
        original_smiles=(row.get('Ligand SMILES') or '').strip()
        smiles=canonical_smiles(original_smiles)
        if not smiles:
            counts['excluded_smiles']+=1;continue
        for endpoint,number,text in labels:
            counts['eligible_measurements']+=1
            key=hashlib.sha256((seq+'\0'+smiles+'\0'+endpoint).encode()).hexdigest()
            seen[key]=seen.get(key,0)+1
            keep=False
            if seen[key]==1:
                counts['distinct_pairs']+=1
                if len(slots)<count:
                    slots.append(key);keep=True
                else:
                    index=rng.randrange(counts['distinct_pairs'])
                    if index<count:
                        del selected[slots[index]]
                        slots[index]=key;keep=True
            elif key in selected:
                keep=label_rng.randrange(seen[key])==0
                counts['repeated_measurements']+=1
            if keep:
                selected[key]=dict(pair_key=key,sequence=seq,smiles=smiles,original_smiles=original_smiles,
                    target_name=row.get('Target Name') or 'Unnamed target',
                    accession=row.get('UniProt (SwissProt) Primary ID of Target Chain 1') or '',
                    ligand_name=row.get('BindingDB MonomerID') or 'Unnamed ligand',
                    endpoint=endpoint,label_nm=number,label_text=text,source_row=row_number,
                    reactant_set_id=row.get('BindingDB Reactant_set_id') or '',
                    article_doi=row.get('Article DOI') or '',pmid=row.get('PMID') or '',
                    publication_date=row.get('Date of publication') or '',
                    curation_date=row.get('Date in BindingDB') or '')
    if progress: progress(dict(counts))
    if len(slots)<count:
        raise ValueError(f'Only {len(slots)} distinct eligible pairs; requested {count}. Reduce the sample size or broaden filters.')
    rng.shuffle(slots)
    return [selected[k] for k in slots],dict(counts)


def snapshot(record,config):
    target_id=hashlib.sha256(record['sequence'].encode()).hexdigest()
    ligand_id=hashlib.sha256(record['smiles'].encode()).hexdigest()
    ref=dict(endpoint=record['endpoint'],value=record['label_text'],unit='nM',source=record['article_doi'])
    return dict(target_id=target_id,target_name=record['target_name'],accession=record['accession'],
        parent_sequence=record['sequence'],sequence=record['sequence'],variant_id='original',variant_name='Original',changes=[],
        ligand_id=ligand_id,ligand_name=record['ligand_name'],smiles=record['smiles'],settings=config['settings'],
        model_version=MODEL_VERSION,model_image='nvcr.io/nim/mit/boltz2:1.9.0',
        alignment='query-only',conversion=CONVERSION,references=[ref],benchmark_label=record)


def ranks(values):
    positions=defaultdict(list)
    for i,value in enumerate(sorted(values),1): positions[value].append(i)
    return [statistics.mean(positions[v]) for v in values]


def correlation(x,y):
    if len(x)<3 or len(set(x))<2 or len(set(y))<2:return None
    return statistics.correlation(x,y)


def statistics_for(points):
    x=[p['label_p'] for p in points];y=[p['predicted_pic50'] for p in points]
    return dict(n=len(points),pearson=correlation(x,y),spearman=correlation(ranks(x),ranks(y)),
        rmse_log10=math.sqrt(statistics.mean((a-b)**2 for a,b in zip(x,y))) if x else None,
        mae_log10=statistics.mean(abs(a-b) for a,b in zip(x,y)) if x else None)


class Benchmarks:
    def __init__(self,store,service):
        self.store=store;self.service=service
        self.lock=threading.Lock()

    def recover(self):
        for row in self.store.all("SELECT run_id FROM benchmarks WHERE state='sampling'"):
            self.prepare_async(row['run_id'])

    def create(self,name,config):
        run_id=uid()
        with self.store.connect() as db:
            db.execute('INSERT INTO runs(id,name,created) VALUES (?,?,?)',(run_id,name,now()))
            db.execute('INSERT INTO benchmarks VALUES (?,?,?,?,?,?)',
                       (run_id,json.dumps(config),'sampling','{}',None,now()))
        self.prepare_async(run_id)
        return self.summary(run_id)

    def prepare_async(self,id):
        threading.Thread(target=self.prepare,args=(id,),daemon=True).start()

    def prepare(self,id):
        with self.lock:
            try:
                campaign=self.store.one('SELECT * FROM benchmarks WHERE run_id=?',(id,))
                if not campaign or campaign['state']!='sampling':return
                config=json.loads(campaign['config'])
                folder=self.store.folder/'data'
                # The bootstrap downloader owns the same cache. Avoid racing its .part file.
                with self.service.bootstrap_lock:
                    prepare_sample(folder)
                archive=folder/'raw'/ARCHIVE
                digest=hashlib.sha256(archive.read_bytes()).hexdigest()
                if hashlib.md5(archive.read_bytes()).hexdigest()!='8c41a1fcf8b828d99070f4bca6bd7a86':
                    raise ValueError('BindingDB release checksum mismatch.')
                def stopped():
                    row=self.store.one('SELECT state FROM benchmarks WHERE run_id=?',(id,))
                    return self.service.stop.is_set() or not row or row['state']=='cancelled'
                def progress(value):
                    with self.store.connect() as db:
                        db.execute('UPDATE benchmarks SET sampling=? WHERE run_id=?',(json.dumps(value),id))
                with zipfile.ZipFile(archive) as z,z.open(MEMBER) as f:
                    reader=csv.DictReader(io.TextIOWrapper(f,encoding='utf-8'),delimiter='\t')
                    selected,counts=sample_rows(reader,config['count'],config['seed'],config['endpoints'],
                        config['min_residues'],config['max_residues'],progress,stopped)
                if stopped():return
                manifest=dict(archive=ARCHIVE,sha256=digest,seed=config['seed'],filters=config,counts=counts,
                    method='Uniform reservoir sampling without replacement over distinct (sequence, canonical SMILES, endpoint) pairs. For each selected pair, one eligible exact measurement is chosen uniformly from its repeated source rows.',
                    limitations='Articles subset, not all BindingDB. Query-only MSA. Assay construct and training overlap unverified. Pooled correlation can reflect between-target effects. Ki/Kd comparisons are proxies, not direct IC50 validation.')
                out=self.store.folder/'benchmarks'/id
                out.mkdir(parents=True,exist_ok=True)
                (out/'manifest.json').write_text(json.dumps(manifest,indent=2))
                with (out/'sample.jsonl').open('w') as f:
                    for record in selected:f.write(json.dumps(record)+'\n')
                with self.store.connect() as db:
                    state=db.execute('SELECT state FROM benchmarks WHERE run_id=?',(id,)).fetchone()['state']
                    if state=='cancelled':return
                    for record in selected:
                        db.execute('INSERT INTO jobs(id,run_id,snapshot,status,created) VALUES (?,?,?,?,?)',
                            (uid(),id,json.dumps(snapshot(record,config)),'queued',now()))
                    db.execute('UPDATE benchmarks SET state=?,sampling=?,error=NULL WHERE run_id=?',
                               ('ready' if state=='sampling' else state,json.dumps(manifest),id))
            except InterruptedError:
                pass
            except Exception as e:
                with self.store.connect() as db:
                    db.execute("UPDATE benchmarks SET state='failed',error=? WHERE run_id=? AND state!='cancelled'",(str(e)[:1000],id))

    def summary(self,id):
        row=self.store.one('SELECT benchmarks.*,runs.name FROM benchmarks JOIN runs ON runs.id=benchmarks.run_id WHERE run_id=?',(id,))
        if not row:return None
        row['config']=json.loads(row['config']);row['sampling']=json.loads(row['sampling'])
        counts={k:0 for k in ['queued','running','succeeded','failed','interrupted','cancelled']}
        for r in self.store.all('SELECT status,COUNT(*) AS n FROM jobs WHERE run_id=? GROUP BY status',(id,)):counts[r['status']]=r['n']
        row['counts']=counts;row['total']=sum(counts.values());row['requested']=row['config']['count']
        state=row['state']
        row['status']=state if state in ('sampling','paused','cancelled','failed') else (
            'running' if counts['running'] else 'queued' if counts['queued'] else 'completed with errors' if counts['failed'] or counts['interrupted'] else 'complete')
        durations=self.store.all("SELECT (julianday(finished)-julianday(started))*86400 AS seconds FROM jobs WHERE run_id=? AND status='succeeded' ORDER BY finished DESC LIMIT 50",(id,))
        avg=statistics.mean(d['seconds'] for d in durations) if durations else None
        row['mean_prediction_seconds']=avg
        row['estimated_remaining_seconds']=avg*(counts['queued']+counts['running']) if avg else None
        row['queue_paused']=self.store.meta('paused',False)
        return row

    def records(self,id):
        # Select only plotting/provenance fields, never thousands of full sequences/CIFs.
        return self.store.all('''SELECT id,status,error,attempt,started,finished,
            json_extract(snapshot,'$.target_name') AS target,
            json_extract(snapshot,'$.target_id') AS target_id,
            json_extract(snapshot,'$.ligand_name') AS ligand,
            json_extract(snapshot,'$.benchmark_label.endpoint') AS endpoint,
            json_extract(snapshot,'$.benchmark_label.label_nm') AS label_nm,
            json_extract(snapshot,'$.benchmark_label.source_row') AS source_row,
            json_extract(snapshot,'$.benchmark_label.reactant_set_id') AS reactant_set_id,
            json_extract(snapshot,'$.benchmark_label.article_doi') AS article_doi,
            json_extract(score,'$.pic50') AS predicted_pic50,
            json_extract(score,'$.ic50_nm') AS predicted_ic50_nm,
            json_extract(score,'$.raw') AS raw_log10_ic50_uM,
            json_extract(score,'$.ligand_iptm') AS ligand_iptm,
            json_extract(score,'$.binder_probability') AS binder_probability
            FROM jobs WHERE run_id=? ORDER BY rowid''',(id,))

    def report(self,id,endpoint='IC50',target=''):
        records=self.records(id)
        points=[]
        for r in records:
            if r['status']!='succeeded' or r['endpoint']!=endpoint or target.lower() not in r['target'].lower():continue
            points.append(dict(r,label_p=9-math.log10(r['label_nm'])))
        grouped=defaultdict(list)
        for p in points:grouped[p['target_id']].append(p)
        within=[]
        for group in grouped.values():
            if len(group)>=5:
                within.append(dict(target=group[0]['target'],**statistics_for(group)))
        stats=statistics_for(points)
        # All successful records contribute to metrics. Bound browser payload only.
        plotted=points if len(points)<=5000 else random.Random(42).sample(points,5000)
        return dict(endpoint=endpoint,metrics=stats,points=plotted,plotted=len(plotted),
            distinct_targets=len(grouped),within_target=within,
            note='Pooled analysis across targets; this is not a held-out benchmark. Labels and model training data may overlap.'
                 +(' Experimental '+endpoint+' is compared with predicted IC50 as a proxy.' if endpoint!='IC50' else ''))

    def action(self,id,action):
        with self.store.connect() as db:
            row=db.execute('SELECT state FROM benchmarks WHERE run_id=?',(id,)).fetchone()
            if not row:raise ValueError('Benchmark not found.')
            if action=='pause':
                if row['state']!='ready':raise ValueError('Only a prepared benchmark can be paused.')
                db.execute("UPDATE benchmarks SET state='paused' WHERE run_id=?",(id,))
            elif action=='resume':
                if row['state']!='paused':raise ValueError('This benchmark is not paused.')
                db.execute("UPDATE benchmarks SET state='ready' WHERE run_id=?",(id,))
            elif action=='cancel':
                db.execute("UPDATE benchmarks SET state='cancelled' WHERE run_id=?",(id,))
                db.execute('UPDATE runs SET cancelled=1 WHERE id=?',(id,))
                db.execute("UPDATE jobs SET status='cancelled',finished=? WHERE run_id=? AND status='queued'",(now(),id))
            elif action=='retry':
                if self.store.meta('paused',False):raise ValueError('Resolve the uncertain backend request and resume the queue first.')
                if row['state']=='failed':
                    db.execute("UPDATE benchmarks SET state='sampling',error=NULL WHERE run_id=?",(id,))
                else:
                    db.execute("UPDATE jobs SET status='queued',error=NULL,score=NULL WHERE run_id=? AND status IN ('failed','interrupted')",(id,))
                    if row['state']=='cancelled':raise ValueError('A cancelled benchmark cannot be resumed. Create a new one.')
            else:raise ValueError('Unknown action.')
        if action=='retry' and row['state']=='failed':self.prepare_async(id)
        return self.summary(id)
