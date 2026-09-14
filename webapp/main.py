"""Local Boltz-2 workspace: static UI, a small JSON API, SQLite, one worker."""
import csv
import io
import json
import os
import sqlite3
from contextlib import asynccontextmanager
from pathlib import Path
from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse, Response
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field
from webapp import domain
from webapp.service import Service
from webapp.benchmark_api import routes as benchmark_routes
from webapp.store import Store, uid, now

ROOT=Path(__file__).resolve().parents[1]


class TargetInput(BaseModel):
    name: str = Field(min_length=1,max_length=100)
    sequence: str = Field(min_length=4,max_length=10000)
    accession: str = Field(default='',max_length=80)


class VariantInput(BaseModel):
    name: str = Field(default='',max_length=100)
    changes: str = Field(min_length=1,max_length=500)


class LigandInput(BaseModel):
    name: str = Field(min_length=1,max_length=100)
    smiles: str = Field(min_length=1,max_length=2000)
    target_id: str | None = None
    endpoint: str = Field(default='',max_length=10)
    value: str = Field(default='',max_length=40)
    source: str = Field(default='',max_length=200)


class CSVInput(BaseModel):
    text: str = Field(max_length=1_000_000)
    target_id: str | None = None


class RunInput(BaseModel):
    name: str = Field(default='Affinity screen',min_length=1,max_length=100)
    target_id: str
    variant_ids: list[str] = Field(min_length=1,max_length=6)
    ligand_ids: list[str] = Field(min_length=1,max_length=50)
    settings: dict = Field(default_factory=dict)


def create_app(folder=None, endpoint=None, background=True, bootstrap=True):
    store=Store(folder or os.getenv('APP_DATA_DIR',str(ROOT/'app-data')))
    service=Service(store,endpoint or os.getenv('BOLTZ_URL','http://boltz2:8000'),bootstrap)

    @asynccontextmanager
    async def lifespan(app):
        if background:
            service.start(worker=os.getenv("RUN_WORKER", "1")=="1")
        yield
        if background:
            service.close()

    app=FastAPI(title='Affinity lab',lifespan=lifespan)
    app.state.store=store
    app.state.service=service

    @app.middleware('http')
    async def local_requests(request: Request, call_next):
        # Browser writes must originate from this same workspace. No CORS or cookies.
        if request.method not in ('GET','HEAD','OPTIONS'):
            origin=request.headers.get('origin')
            if origin and origin != str(request.base_url).rstrip('/'):
                return JSONResponse({'detail':'Cross-origin writes are not allowed.'},status_code=403)
            length=request.headers.get('content-length','0')
            if length.isdigit() and int(length)>1_100_000:
                return JSONResponse({'detail':'Request exceeds 1 MB.'},status_code=413)
        response=await call_next(request)
        response.headers['X-Content-Type-Options']='nosniff'
        response.headers['Referrer-Policy']='same-origin'
        return response

    @app.exception_handler(ValueError)
    async def validation_error(request, exc):
        return JSONResponse({'detail':str(exc)},status_code=422)

    def get_target(id):
        target=store.one('SELECT * FROM targets WHERE id=?',(id,))
        if not target: raise HTTPException(404,'Target not found.')
        return target

    def get_run(id):
        result=store.run(id)
        if result is None: raise HTTPException(404,'Run not found.')
        return result

    @app.get('/')
    def index():
        return FileResponse(ROOT/'webapp/static/index.html')

    @app.get('/api/status')
    def status():
        return dict(ready=service.ready,message=service.health_message,paused=store.meta('paused',False),
                    bootstrap=store.meta('bootstrap',{'status':'idle'}),worker=store.meta('worker',{}),model='Boltz-2',version=domain.MODEL_VERSION,
                    active=store.one("SELECT id,run_id FROM jobs WHERE status='running'"))

    @app.post('/api/bootstrap')
    def retry_bootstrap():
        service.bootstrap_async()
        return {'ok':True}

    @app.post('/api/queue/resume')
    def resume(data: dict):
        if data.get('backend_restarted') is not True:
            raise HTTPException(409,'Restart the Boltz-2 service before resuming uncertain requests.')
        if not service.ready:
            raise HTTPException(409,'Wait for the restarted model to become ready.')
        store.set_meta('paused',False)
        return {'ok':True}

    @app.get('/api/targets')
    def targets():
        return store.all('SELECT * FROM targets ORDER BY name')

    @app.post('/api/targets')
    def add_target(data: TargetInput):
        seq=domain.sequence(data.sequence)
        if not data.name.strip(): raise ValueError('Give the protein a name.')
        id=uid()
        with store.connect() as db:
            db.execute('INSERT INTO targets VALUES (?,?,?,?)',(id,data.name.strip(),seq,data.accession.strip()))
        return get_target(id)

    @app.get('/api/targets/{id}/variants')
    def variants(id: str):
        get_target(id)
        rows=store.all('SELECT * FROM variants WHERE target_id=? ORDER BY rowid',(id,))
        for row in rows: row['changes']=json.loads(row['changes'])
        return rows

    @app.post('/api/targets/{id}/preview')
    def preview(id: str,data: VariantInput):
        target=get_target(id)
        seq,changes=domain.mutate(target['sequence'],data.changes)
        return dict(sequence=seq,changes=changes)

    @app.post('/api/targets/{id}/variants')
    def add_variant(id: str,data: VariantInput):
        result=preview(id,data)
        vid=uid()
        name=data.name.strip() or ', '.join(c['label'] for c in result['changes'])
        with store.connect() as db:
            db.execute('INSERT INTO variants VALUES (?,?,?,?,?)',
                (vid,id,name,result['sequence'],json.dumps(result['changes'])))
        return dict(id=vid,target_id=id,name=name,**result)

    @app.get('/api/ligands')
    def ligands():
        rows=store.all('SELECT id,name,smiles,mw FROM ligands ORDER BY rowid')
        refs=store.all('SELECT * FROM measurements')
        for row in rows:
            row['references']=[r for r in refs if r['ligand_id']==row['id']]
        return rows

    @app.get('/api/ligands/{id}/image')
    def ligand_image(id: str):
        row=store.one('SELECT svg FROM ligands WHERE id=?',(id,))
        if not row: raise HTTPException(404,'Ligand not found.')
        return Response(row['svg'],media_type='image/svg+xml',headers={'Cache-Control':'public,max-age=86400'})

    def insert_ligand(data):
        smiles,svg,mw=domain.molecule(data.smiles)
        if not data.name.strip(): raise ValueError('Give the ligand a name.')
        if data.target_id: get_target(data.target_id)
        if data.value and not data.target_id: raise ValueError('Select the target for this reference measurement.')
        if data.value and data.endpoint not in ('IC50','Ki','Kd','EC50'):
            raise ValueError('Reference endpoint must be IC50, Ki, Kd, or EC50.')
        if data.value:
            import re
            if not re.fullmatch(r'[<>~]?\s*\d*\.?\d+(?:[eE][+-]?\d+)?',data.value.strip()):
                raise ValueError('Reference value must be a positive number, optionally with < or >, in nM.')
            number=float(re.sub(r'^[<>~]\s*','',data.value.strip()))
            if not 0<number<1e100: raise ValueError('Reference value must be positive and finite.')
        with store.connect() as db:
            existing=db.execute('SELECT id FROM ligands WHERE smiles=?',(smiles,)).fetchone()
            id=existing['id'] if existing else uid()
            if not existing:
                db.execute('INSERT INTO ligands VALUES (?,?,?,?,?)',(id,data.name.strip(),smiles,svg,mw))
            if data.target_id:
                db.execute('INSERT OR IGNORE INTO measurements VALUES (?,?,?,?,?,?)',
                    (data.target_id,id,data.endpoint,data.value.strip(),'nM' if data.value else '',data.source))
        return dict(id=id,name=data.name.strip(),smiles=smiles,mw=mw,existing=bool(existing))

    @app.post('/api/ligands')
    def add_ligand(data: LigandInput):
        return insert_ligand(data)

    @app.post('/api/ligands/import')
    def import_ligands(data: CSVInput):
        reader=csv.DictReader(io.StringIO(data.text))
        if not {'name','smiles'}.issubset(reader.fieldnames or []):
            raise ValueError('CSV needs name and smiles columns.')
        rows=list(reader)
        if len(rows)>500: raise ValueError('Import at most 500 rows at a time.')
        added,errors=[],[]
        for number,row in enumerate(rows,2):
            try:
                if row.get('unit') and row['unit']!='nM':
                    raise ValueError('Reference measurements must be in nM.')
                item=LigandInput(name=row['name'],smiles=row['smiles'],target_id=data.target_id,
                    endpoint=row.get('endpoint') or '',value=row.get('value') or '',source=row.get('source') or '')
                added.append(insert_ligand(item))
            except (ValueError,TypeError) as e:
                errors.append(dict(row=number,message=str(e)[:250]))
        return dict(imported=added,errors=errors)

    @app.post('/api/runs')
    def create_run(data: RunInput):
        if not service.ready: raise HTTPException(409,'Boltz-2 is not ready yet.')
        if store.meta('paused',False): raise HTTPException(409,'Resolve the interrupted request before starting another run.')
        target=get_target(data.target_id)
        chosen=list(dict.fromkeys(data.variant_ids))
        ligand_ids=list(dict.fromkeys(data.ligand_ids))
        if len(chosen)*len(ligand_ids)>100: raise ValueError('Limit a run to 100 predictions.')
        config=domain.settings(data.settings)
        selected_variants=[]
        for vid in chosen:
            if vid=='original':
                selected_variants.append(dict(id='original',name='Original',sequence=target['sequence'],changes=[]))
            else:
                variant=store.one('SELECT * FROM variants WHERE id=? AND target_id=?',(vid,target['id']))
                if not variant: raise ValueError('A selected variant does not belong to this target.')
                variant['changes']=json.loads(variant['changes'])
                selected_variants.append(variant)
        selected_ligands=[]
        for lid in ligand_ids:
            ligand=store.one('SELECT id,name,smiles FROM ligands WHERE id=?',(lid,))
            if not ligand: raise ValueError('A selected ligand no longer exists.')
            selected_ligands.append(ligand)
        id=uid()
        with store.connect() as db:
            db.execute('INSERT INTO runs(id,name,created) VALUES (?,?,?)',(id,data.name.strip() or 'Affinity screen',now()))
            for variant in selected_variants:
                for ligand in selected_ligands:
                    refs=[dict(r) for r in db.execute('SELECT * FROM measurements WHERE target_id=? AND ligand_id=?',
                                                    (target['id'],ligand['id']))] if variant['id']=='original' else []
                    snapshot=dict(target_id=target['id'],target_name=target['name'],accession=target['accession'],
                        parent_sequence=target['sequence'],sequence=variant['sequence'],variant_id=variant['id'],
                        variant_name=variant['name'],changes=variant['changes'],ligand_id=ligand['id'],ligand_name=ligand['name'],
                        smiles=ligand['smiles'],settings=config,model_version=domain.MODEL_VERSION,
                        model_image='nvcr.io/nim/mit/boltz2:1.9.0',
                        model_digest='sha256:a2eb69458fd9505ce9fdaebb16e206da909666ec8d94a297d81e7b57a2f6c090',
                        alignment='query-only',conversion=domain.CONVERSION,references=refs)
                    db.execute('INSERT INTO jobs(id,run_id,snapshot,status,created) VALUES (?,?,?,?,?)',
                               (uid(),id,json.dumps(snapshot),'queued',now()))
        return get_run(id)

    @app.get('/api/runs')
    def runs():
        result=[]
        for row in store.all('SELECT id FROM runs WHERE id NOT IN (SELECT run_id FROM benchmarks) ORDER BY created DESC LIMIT 100'):
            run=get_run(row['id'])
            run['target_name']=run['jobs'][0]['snapshot']['target_name'] if run['jobs'] else ''
            del run['jobs']
            result.append(run)
        return result

    @app.get('/api/runs/{id}')
    def read_run(id: str):
        return get_run(id)

    @app.post('/api/runs/{id}/cancel')
    def cancel(id: str):
        get_run(id)
        with store.connect() as db:
            db.execute('UPDATE runs SET cancelled=1 WHERE id=?',(id,))
            db.execute("UPDATE jobs SET status='cancelled',finished=? WHERE run_id=? AND status='queued'",(now(),id))
        return get_run(id)

    @app.post('/api/runs/{id}/retry')
    def retry(id: str):
        get_run(id)
        if store.meta('paused',False): raise HTTPException(409,'Restart the backend and resume the queue first.')
        with store.connect() as db:
            db.execute('UPDATE runs SET cancelled=0 WHERE id=?',(id,))
            db.execute("UPDATE jobs SET status='queued',error=NULL,score=NULL WHERE run_id=? AND status IN ('failed','interrupted')",(id,))
        return get_run(id)

    @app.get('/api/jobs/{id}/attempts')
    def attempts(id: str):
        return store.all('SELECT * FROM attempts WHERE job_id=? ORDER BY number',(id,))

    @app.get('/api/jobs/{id}/artifacts/{kind}')
    def artifact(id: str,kind: str):
        names={'structure':('complex.cif','chemical/x-mmcif'), 'response':('response.json','application/json'),
               'request':('request.json','application/json'),'scores':('scores.json','application/json'),
               'error':('error.txt','text/plain')}
        if kind not in names: raise HTTPException(404,'Artifact type not found.')
        job=store.one('SELECT attempt FROM jobs WHERE id=?',(id,))
        if not job: raise HTTPException(404,'Job not found.')
        name,mime=names[kind]
        path=store.folder/'artifacts'/id/str(job['attempt'])/name
        if not path.is_file(): raise HTTPException(404,'This artifact is not available.')
        return FileResponse(path,media_type=mime,filename=f'{id}-{name}')

    @app.get('/api/runs/{id}/export')
    def export(id: str):
        run=get_run(id)
        out=io.StringIO()
        fields=['target','variant','mutations','ligand','smiles','status','rank','predicted_ic50_equivalent_nm','predicted_pic50','delta_pic50_vs_original','raw_log10_ic50_uM','binder_probability','ligand_iptm','model_version','alignment','ground_truth_endpoint','ground_truth_value','ground_truth_unit','ground_truth_source','reference','error']
        writer=csv.DictWriter(out,fieldnames=fields)
        writer.writeheader()
        baseline={j['snapshot']['ligand_id']:j['score']['pic50'] for j in run['jobs'] if j['status']=='succeeded' and j['snapshot']['variant_id']=='original'}
        for j in sorted(run['jobs'],key=lambda j:(j['snapshot']['variant_name'],j['rank'] if j['rank'] is not None else 1e9)):
            s=j['snapshot']; v=j['score'] or {}; base=baseline.get(s['ligand_id'])
            row=dict(target=s['target_name'],variant=s['variant_name'],mutations=','.join(c['label'] for c in s['changes']),
                ligand=s['ligand_name'],smiles=s['smiles'],status=j['status'],rank=j['rank'],predicted_ic50_equivalent_nm=v.get('ic50_nm'),
                predicted_pic50=v.get('pic50'),delta_pic50_vs_original=v['pic50']-base if base is not None and v else None,
                raw_log10_ic50_uM=v.get('raw'),binder_probability=v.get('binder_probability'),ligand_iptm=v.get('ligand_iptm'),
                model_version=s['model_version'],alignment=s['alignment'],
                ground_truth_endpoint='; '.join(r['endpoint'] for r in s['references'] if r['value']),
                ground_truth_value='; '.join(r['value'] for r in s['references'] if r['value']),
                ground_truth_unit='; '.join(r['unit'] for r in s['references'] if r['value']),
                ground_truth_source='; '.join(r['source'] for r in s['references'] if r['value']),
                reference=json.dumps(s['references']),error=j['error'])
            # Prevent imported labels from becoming formulas in spreadsheet software.
            writer.writerow({k:("'"+value if isinstance(value,str) and value.startswith(('=','+','-','@','\t','\r')) else value) for k,value in row.items()})
        return Response(out.getvalue(),media_type='text/csv',headers={'Content-Disposition':f'attachment; filename="ranking-{id}.csv"'})

    benchmark_routes(app,store,service)

    app.mount('/static',StaticFiles(directory=ROOT/'webapp/static'),name='static')
    app.mount('/vendor',StaticFiles(directory=ROOT/'assets/vendor'),name='vendor')
    return app


app=create_app()
