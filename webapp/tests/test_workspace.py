import io
import csv
import json
import urllib.error
from pathlib import Path

import pytest
from fastapi.testclient import TestClient
from webapp.domain import mutate, scores, sequence
from webapp.main import create_app


@pytest.fixture
def workspace(tmp_path):
    app=create_app(tmp_path,background=False,bootstrap=False)
    app.state.service.ready=True
    with TestClient(app) as client:
        target=client.post('/api/targets',json={'name':'Test protein','sequence':'>protein\nACDEFGHIK'}).json()
        first=client.post('/api/ligands',json={'name':'Ligand A','smiles':'CCO','target_id':target['id'],'endpoint':'IC50','value':'>10'}).json()
        second=client.post('/api/ligands',json={'name':'Ligand B','smiles':'CCN'}).json()
        yield app,client,target,[first,second]


def create_run(client,target,ligands,variant_ids=None):
    return client.post('/api/runs',json={'name':'Test screen','target_id':target['id'],
        'variant_ids':variant_ids or ['original'],'ligand_ids':[x['id'] for x in ligands]})


def answer(raw=-1):
    return {'affinities':{'B':{'affinity_pred_value':[raw],'affinity_pic50':[99],
        'affinity_probability_binary':[0.7]}},'ligand_iptm_scores':[0.9],
        'structures':[{'structure':'data_test\n#\n'}]}


def mock_response(monkeypatch, raw=-1):
    def fake(req,timeout):
        request=json.loads(req.data)
        assert len(request['polymers'])==1
        assert len(request['ligands'])==1
        assert request['ligands'][0]['predict_affinity'] is True
        return io.BytesIO(json.dumps(answer(raw)).encode())
    monkeypatch.setattr('urllib.request.urlopen',fake)


def test_mutations_are_relative_to_original_and_preserve_length():
    result,changes=mutate('ACDEFGHIK','C2Y, I8L')
    assert result=='AYDEFGHLK'
    assert len(result)==9
    assert [c['label'] for c in changes]==['C2Y','I8L']
    for bad in ['C0Y','C10Y','T2Y','C2C','C2Y C2F','C2Z','del2']:
        with pytest.raises(ValueError): mutate('ACDEFGHIK',bad)
    with pytest.raises(ValueError): sequence('>one\nACDE\n>two\nACDE')


def test_affinity_conversion_ignores_mislabeled_nim_field():
    result=scores(answer(-1.53125))
    assert result['pic50']==7.53125
    assert result['ic50_nm']==pytest.approx(29.427271762)
    for bad in [float('nan'),float('inf'),True,'bad',1000,-1000]:
        with pytest.raises(ValueError): scores(answer(bad))
    with pytest.raises(ValueError): scores(answer(), 'unknown')


def test_run_cross_product_snapshots_and_refs(workspace):
    app,c,t,ligands=workspace
    v=c.post(f"/api/targets/{t['id']}/variants",json={'changes':'C2Y'}).json()
    r=create_run(c,t,ligands,['original',v['id']]).json()
    assert r['total']==4
    assert [j['snapshot']['sequence'] for j in r['jobs']]==[t['sequence'],t['sequence'],'AYDEFGHIK','AYDEFGHIK']
    assert r['jobs'][0]['snapshot']['references'][0]['value']=='>10'
    assert r['jobs'][2]['snapshot']['references']==[]
    exported=list(csv.DictReader(io.StringIO(c.get('/api/runs/'+r['id']+'/export').text)))
    original=[row for row in exported if row['variant']=='Original' and row['ligand']=='Ligand A'][0]
    assert original['ground_truth_value']=='>10' and original['ground_truth_endpoint']=='IC50'
    assert all(not row['ground_truth_value'] for row in exported if row['variant']!='Original')
    with app.state.store.connect() as db:
        db.execute('UPDATE targets SET sequence=? WHERE id=?',('YYYY',t['id']))
    assert c.get('/api/runs/'+r['id']).json()['jobs'][0]['snapshot']['sequence']==t['sequence']
    app.state.service.ready=False
    assert create_run(c,t,ligands).status_code==409


def test_serial_claim_cancel_and_tied_ranking(workspace,monkeypatch):
    app,c,t,ligands=workspace
    run=create_run(c,t,ligands).json()
    mock_response(monkeypatch)
    first=app.state.store.claim()
    app.state.service.predict(first)
    second=app.state.store.claim()
    app.state.service.predict(second)
    result=c.get('/api/runs/'+run['id']).json()
    assert result['status']=='complete'
    assert [j['rank'] for j in result['jobs']]==[1.5,1.5]
    assert '100.0' in c.get('/api/runs/'+run['id']+'/export').text
    assert c.get('/api/jobs/'+first['id']+'/artifacts/structure').status_code==200
    assert c.get('/api/jobs/'+first['id']+'/artifacts/../../workspace.sqlite3').status_code==404
    next_run=create_run(c,t,ligands).json()
    active=app.state.store.claim()
    c.post('/api/runs/'+next_run['id']+'/cancel')
    result=c.get('/api/runs/'+next_run['id']).json()
    assert result['counts']['running']==1 and result['counts']['cancelled']==1
    app.state.service.predict(active)
    assert app.state.store.claim() is None


def test_failed_job_keeps_success_and_retry_attempts(workspace,monkeypatch):
    app,c,t,ligands=workspace
    run=create_run(c,t,ligands).json()
    job=app.state.store.claim()
    def fail(*args,**kwargs):
        raise urllib.error.HTTPError('http://test',422,'bad',{},io.BytesIO(b'invalid input'))
    monkeypatch.setattr('urllib.request.urlopen',fail)
    app.state.service.predict(job)
    mock_response(monkeypatch)
    second=app.state.store.claim()
    app.state.service.predict(second)
    result=c.get('/api/runs/'+run['id']).json()
    assert result['status']=='completed with errors'
    assert result['counts']['failed']==1 and result['counts']['succeeded']==1
    c.post('/api/runs/'+run['id']+'/retry')
    retried=app.state.store.claim()
    assert retried['id']==job['id'] and retried['attempt']==2
    app.state.service.predict(retried)
    assert c.get('/api/runs/'+run['id']).json()['status']=='complete'
    assert len(c.get('/api/jobs/'+job['id']+'/attempts').json())==2
    assert (app.state.store.folder/'artifacts'/job['id']/'1'/'error.txt').exists()


def test_restart_pauses_uncertain_jobs_and_keeps_queue(workspace):
    app,c,t,ligands=workspace
    run=create_run(c,t,ligands).json()
    app.state.store.claim()
    app.state.store.recover()
    result=c.get('/api/runs/'+run['id']).json()
    assert result['counts']['interrupted']==1 and result['counts']['queued']==1
    assert app.state.store.claim() is None
    assert c.post('/api/queue/resume',json={}).status_code==409
    assert c.post('/api/queue/resume',json={'backend_restarted':True}).status_code==200
    assert app.state.store.claim() is not None


def test_library_validation_and_csv_errors(workspace):
    app,c,t,ligands=workspace
    assert c.post('/api/ligands',json={'name':'bad','smiles':'not-a-molecule'}).status_code==422
    assert c.post('/api/ligands',json={'name':'salt','smiles':'CC.[Na+]'}).status_code==422
    assert c.post('/api/targets',json={'name':'bad','sequence':'QQQ*X'}).status_code==422
    result=c.post('/api/ligands/import',json={'text':'name,smiles\nvalid,CCC\nbad,invalid','target_id':t['id']}).json()
    assert len(result['imported'])==1 and result['errors'][0]['row']==3
    assert c.post('/api/runs',json={'target_id':t['id'],'variant_ids':['invalid'],'ligand_ids':[ligands[0]['id']]}).status_code==422
    assert c.post('/api/targets',headers={'Origin':'https://outside.example'},json={'name':'bad','sequence':'ACDE'}).status_code==403
