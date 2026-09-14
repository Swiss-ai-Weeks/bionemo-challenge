"""HTTP routes for long-running BindingDB evaluation campaigns."""
import csv
import io
import json
from fastapi import HTTPException, Query
from fastapi.responses import Response, FileResponse
from pydantic import BaseModel, Field
from webapp.domain import settings


class BenchmarkInput(BaseModel):
    name: str = Field(default='BindingDB affinity benchmark',min_length=1,max_length=100)
    count: int = Field(default=2000,ge=1,le=50000)
    seed: int = Field(default=20260914,ge=0,le=4294967295)
    endpoints: list[str] = Field(default_factory=lambda:['IC50'],min_length=1,max_length=3)
    min_residues: int = Field(default=50,ge=4,le=2048)
    max_residues: int = Field(default=1000,ge=4,le=2048)
    settings: dict = Field(default_factory=dict)


def routes(app,store,service):
    bench=service.benchmarks
    def require(id):
        row=bench.summary(id)
        if not row:raise HTTPException(404,'Benchmark not found.')
        return row

    @app.post('/api/benchmarks')
    def create(data: BenchmarkInput):
        if not set(data.endpoints).issubset({'IC50','Ki','Kd'}):raise ValueError('Choose IC50, Ki, or Kd.')
        if data.min_residues>data.max_residues:raise ValueError('Minimum protein length exceeds maximum.')
        config=data.model_dump(exclude={'name'})
        config['settings']=settings(data.settings)
        config['endpoints']=list(dict.fromkeys(data.endpoints))
        return bench.create(data.name.strip() or 'BindingDB affinity benchmark',config)

    @app.get('/api/benchmarks')
    def list_benchmarks():
        return [bench.summary(r['run_id']) for r in store.all('SELECT run_id FROM benchmarks ORDER BY created DESC LIMIT 100')]

    @app.get('/api/benchmarks/{id}')
    def summary(id: str):return require(id)

    @app.post('/api/benchmarks/{id}/{action}')
    def action(id: str,action: str):
        require(id)
        return bench.action(id,action)

    @app.get('/api/benchmarks/{id}/report')
    def report(id: str,endpoint: str='IC50',target: str=Query(default='',max_length=200)):
        require(id)
        if endpoint not in ('IC50','Ki','Kd'):raise ValueError('Invalid endpoint.')
        return bench.report(id,endpoint,target)

    @app.get('/api/benchmarks/{id}/records')
    def records(id: str,offset: int=Query(default=0,ge=0),limit: int=Query(default=50,ge=1,le=100),
                endpoint: str='',status: str=''):
        require(id)
        rows=bench.records(id)
        rows=[r for r in rows if (not endpoint or r['endpoint']==endpoint) and (not status or r['status']==status)]
        return dict(total=len(rows),offset=offset,records=rows[offset:offset+limit])

    @app.get('/api/benchmarks/{id}/export')
    def export(id: str):
        require(id)
        rows=bench.records(id)
        if not rows:raise HTTPException(409,'Sampling has not finished.')
        output=io.StringIO()
        writer=csv.DictWriter(output,fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow({k:("'"+v if isinstance(v,str) and v.startswith(('=','+','-','@','\t','\r')) else v) for k,v in row.items()})
        return Response(output.getvalue(),media_type='text/csv',headers={'Content-Disposition':f'attachment; filename="benchmark-{id}.csv"'})

    @app.get('/api/benchmarks/{id}/manifest')
    def manifest(id: str):
        require(id)
        path=store.folder/'benchmarks'/id/'manifest.json'
        if not path.is_file():raise HTTPException(409,'Sampling has not finished.')
        return FileResponse(path,media_type='application/json',filename=f'benchmark-{id}-manifest.json')

    @app.get('/api/benchmarks/{id}/sample')
    def sample(id: str):
        require(id)
        path=store.folder/'benchmarks'/id/'sample.jsonl'
        if not path.is_file():raise HTTPException(409,'Sampling has not finished.')
        return FileResponse(path,media_type='application/x-ndjson',filename=f'benchmark-{id}-sample.jsonl')

    @app.get('/api/benchmarks/{id}/plot.svg')
    def plot(id: str,endpoint: str='IC50',target: str=Query(default='',max_length=200)):
        from matplotlib.figure import Figure
        require(id)
        if endpoint not in ('IC50','Ki','Kd'):raise ValueError('Invalid endpoint.')
        report=bench.report(id,endpoint,target)
        fig=Figure(figsize=(7.4,6),layout='constrained')
        ax=fig.subplots()
        points=report['points']
        if points:
            x=[p['label_p'] for p in points];y=[p['predicted_pic50'] for p in points]
            low=min(x+y)-.5;high=max(x+y)+.5
            ax.scatter(x,y,s=12,c='#087f69',alpha=.45,edgecolors='none')
            ax.plot([low,high],[low,high],color='#a1aaa6',linestyle='--',linewidth=1)
            ax.set_xlim(low,high);ax.set_ylim(low,high)
        else:
            ax.text(.5,.5,'Waiting for successful predictions',ha='center',transform=ax.transAxes)
        ax.set_xlabel('BindingDB experimental p'+endpoint+' (−log10 M)')
        ax.set_ylabel('Boltz-2 predicted pIC50')
        ax.set_title(f"BindingDB vs Boltz-2 · {endpoint} · n = {report['metrics']['n']}")
        ax.grid(alpha=.15);ax.set_axisbelow(True)
        fig.text(.02,.005,'Query-only alignment. Pooled exploratory analysis; training overlap and assay constructs unverified.',fontsize=7,color='#65736b')
        output=io.StringIO();fig.savefig(output,format='svg')
        return Response(output.getvalue(),media_type='image/svg+xml',headers={'Content-Disposition':f'attachment; filename="benchmark-{id}-{endpoint}.svg"'})
