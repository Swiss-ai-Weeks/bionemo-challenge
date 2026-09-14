"""Build an offline HTML viewer for a completed DiffDock test (latest by default)."""
import argparse
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('directory', nargs='?', type=Path)
a=p.parse_args()
folder=a.directory or sorted((ROOT/'results/diffdock').glob('*/summary.json'))[-1].parent
summary=json.loads((folder/'summary.json').read_text())
run=summary['runs'][-1]
poses=[{'sdf':f.read_text(),'confidence':score} for f,score in zip(sorted((folder/f"run_{run['run']}").glob('pose_*.sdf')),run['confidence_descending'])]
data=json.dumps({'protein':(folder/'target.pdb').read_text(),'poses':poses}).replace('<','\\u003c')
library=(ROOT/'assets/vendor/3Dmol-min.js').read_text()
if '</script' in library.lower():raise ValueError('Unexpected script tag in library')
html='''<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>DiffDock · Src + 4521</title>
<style>body{font:16px system-ui;margin:24px;color:#172033;background:#f6f8fa}#viewer{height:70vh;position:relative}button,select{padding:8px;margin:8px}p{max-width:900px}</style>
<script>LIBRARY</script></head><body><h1>DiffDock · Src + ligand 4521</h1>
<p>Fixed receptor: Boltz-2 predicted Src. Orange sticks: DiffDock ligand pose. Confidence ranks poses; it is not affinity or experimental validation.</p>
<label>Pose <select id="pose"></select></label><button id="whole">Whole complex</button><button id="ligand">Focus ligand</button><div id="viewer"></div>
<script>const data=DATA;const select=document.getElementById('pose');
data.poses.forEach((p,i)=>{const o=document.createElement('option');o.value=i;o.textContent=`${i+1} · confidence ${p.confidence.toFixed(3)}`;select.append(o);});
const viewer=$3Dmol.createViewer(document.getElementById('viewer'),{backgroundColor:'white'});
function show(){viewer.removeAllModels();viewer.addModel(data.protein,'pdb');viewer.addModel(data.poses[Number(select.value)].sdf,'sdf');viewer.setStyle({model:0},{cartoon:{color:'#4477aa'}});viewer.setStyle({model:1},{stick:{colorscheme:'orangeCarbon',radius:0.22}});viewer.zoomTo();viewer.render();}
select.onchange=show;document.getElementById('whole').onclick=()=>{viewer.zoomTo();viewer.render();};document.getElementById('ligand').onclick=()=>{viewer.zoomTo({model:1});viewer.render();};window.addEventListener('resize',()=>{viewer.resize();viewer.render();});show();
</script></body></html>'''.replace('DATA',data).replace('LIBRARY',library)
(folder/'viewer.html').write_text(html)
print(folder/'viewer.html')
