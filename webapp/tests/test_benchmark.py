import csv
import io
import json
import pytest
from fastapi.testclient import TestClient
from webapp.benchmark import sample_rows, statistics_for, snapshot
from webapp.domain import DEFAULT_SETTINGS
from webapp.main import create_app
from webapp.store import now,uid


def row(i, value=None, endpoint='IC50', seq='ACDEFGHIK'):
    return {'BindingDB Target Chain Sequence 1':seq,'Ligand SMILES':'C'*(i+1),
        'Number of Protein Chains in Target (>1 implies a multichain complex)':'1',
        endpoint+' (nM)':str(value if value is not None else i+1),
        'Target Name':'Protein A','BindingDB MonomerID':str(i),'BindingDB Reactant_set_id':str(i*10)}


def test_sampler_exact_labels_deduplicates_and_is_reproducible():
    source=[row(i) for i in range(10)]
    source += [row(0,10),row(0,100),row(3,'>500'),row(3,'~5'),row(1,0),row(1,'nan')]
    source += [row(1,10,seq='XXXX')]
    one,counts=sample_rows(iter(source),5,42,['IC50'],4,100)
    two,_=sample_rows(iter(source),5,42,['IC50'],4,100)
    assert one==two
    assert len({r['pair_key'] for r in one})==5
    assert counts['distinct_pairs']==10
    assert counts['excluded_censored_or_invalid_labels']==4
    assert all(r['label_nm']>0 and r['endpoint']=='IC50' for r in one)
    with pytest.raises(ValueError):sample_rows(iter(source),11,42,['IC50'],4,100)
    # Selection is uniform over pairs, not heavily repeated source rows.
    repeated=[row(0)]*1000+[row(i) for i in range(1,10)]
    hits=0
    for seed in range(300):
        selected,_=sample_rows(iter(repeated),1,seed,['IC50'],4,100)
        hits+=selected[0]['ligand_name']=='0'
    assert 10<hits<60


def test_statistics_and_endpoint_separation():
    pts=[dict(label_p=x,predicted_pic50=x) for x in [4.,5.,6.]]
    stats=statistics_for(pts)
    assert stats['pearson']==pytest.approx(1)
    assert stats['spearman']==pytest.approx(1)
    assert stats['rmse_log10']==0
    assert statistics_for(pts[:1])['pearson'] is None
    selected,counts=sample_rows(iter([row(0,1),row(0,10,'Ki')]),2,1,['IC50','Ki'],4,100)
    assert {r['endpoint'] for r in selected}=={'IC50','Ki'}


@pytest.fixture
def benchmark(tmp_path,monkeypatch):
    app=create_app(tmp_path,background=False,bootstrap=False)
    app.state.service.ready=True
    monkeypatch.setattr(app.state.service.benchmarks,'prepare_async',lambda id:None)
    with TestClient(app) as c:
        response=c.post('/api/benchmarks',json={'count':3,'min_residues':4,'max_residues':100})
        assert response.status_code==200
        id=response.json()['run_id']
        config=response.json()['config']
        selected,_=sample_rows(iter([row(0,1),row(1,10),row(2,100)]),3,1,['IC50'],4,100)
        with app.state.store.connect() as db:
            for r in selected:
                db.execute('INSERT INTO jobs(id,run_id,snapshot,status,created) VALUES (?,?,?,?,?)',
                    (uid(),id,json.dumps(snapshot(r,config)),'queued',now()))
            db.execute("UPDATE benchmarks SET state='ready' WHERE run_id=?",(id,))
        yield app,c,id


def test_campaign_priority_pause_resume_report_and_export(benchmark):
    app,c,id=benchmark
    assert c.get('/api/runs').json()==[]  # separate from hand-built screens
    t=c.post('/api/targets',json={'name':'Interactive','sequence':'ACDE'}).json()
    l=c.post('/api/ligands',json={'name':'Interactive ligand','smiles':'CCO'}).json()
    r=c.post('/api/runs',json={'target_id':t['id'],'ligand_ids':[l['id']],'variant_ids':['original']}).json()
    claimed=app.state.store.claim()
    assert claimed['run_id']==r['id']
    assert app.state.store.claim() is None  # cannot overlap an active job
    app.state.store.finish(claimed,'failed',error='fixture')
    c.post('/api/benchmarks/'+id+'/pause')
    assert app.state.store.claim() is None
    c.post('/api/benchmarks/'+id+'/resume')
    while job:=app.state.store.claim():
        import math
        nm=json.loads(job['snapshot'])['benchmark_label']['label_nm']
        app.state.store.finish(job,'succeeded',score={'pic50':9-math.log10(nm),'ic50_nm':nm,'raw':math.log10(nm/1000),'binder_probability':.8,'ligand_iptm':.9})
    result=c.get('/api/benchmarks/'+id+'/report').json()
    assert result['metrics']['n']==3
    assert result['metrics']['rmse_log10']==0
    assert c.get('/api/benchmarks/'+id+'/report?endpoint=Ki').json()['metrics']['n']==0
    records=c.get('/api/benchmarks/'+id+'/records?limit=2&offset=1').json()
    assert records['total']==3 and len(records['records'])==2
    rows=list(csv.DictReader(io.StringIO(c.get('/api/benchmarks/'+id+'/export').text)))
    assert len(rows)==3 and all(r['status']=='succeeded' for r in rows)
    svg=c.get('/api/benchmarks/'+id+'/plot.svg')
    assert svg.status_code==200 and '<svg' in svg.text


def test_web_restart_does_not_interrupt_worker_job(benchmark,monkeypatch):
    app,c,id=benchmark
    active=app.state.store.claim()
    monkeypatch.setenv('RUN_WORKER','0')
    restarted=create_app(app.state.store.folder,endpoint='http://127.0.0.1:1',background=True,bootstrap=False)
    with TestClient(restarted):
        record=restarted.state.store.one('SELECT status FROM jobs WHERE id=?',(active['id'],))
        assert record['status']=='running'
        assert not restarted.state.store.meta('paused',False)
        assert restarted.state.store.one('SELECT COUNT(*) AS n FROM jobs WHERE run_id=?',(id,))['n']==3
