"""Create an HTML complex viewer with embedded predictions (no structure upload)."""
import json
from pathlib import Path
root = Path(__file__).resolve().parents[1]
models = {}
for model in ['openfold3', 'boltz2']:
    p = root / 'results' / model / 'complex_0.cif'
    if p.exists():
        models[model] = {'cif': p.read_text(), 'scores': json.loads((p.parent / 'scores.json').read_text())}
        summary = p.parent / 'affinity_summary.json'
        if summary.exists():
            models[model]['affinity'] = json.loads(summary.read_text())
if not models:
    raise SystemExit('Run a prediction first')
comparison = None
if len(models) == 2:
    from compare_structures import compare
    comparison = compare()
    models['openfold3']['cif'] = (root / 'results/comparison/openfold3_aligned.cif').read_text()
data = json.dumps(models).replace('<', '\\u003c')
library = (root / 'assets' / 'vendor' / '3Dmol-min.js').read_text()
if '</script' in library.lower():
    raise ValueError('Bundled library contains an HTML script closing tag')
html = '''<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Src–4521 complexes</title>
<script>LIBRARYDATA</script>
<style>
body{font:16px system-ui;margin:24px;background:#f6f8fa;color:#172033}
.comparison{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:24px}
.panel{min-width:0;background:white;border:1px solid #d8dee6;border-radius:12px;padding:20px}
h2{margin-top:0}.viewer{height:55vh;min-height:320px;position:relative}
.metrics{min-height:4.5em;line-height:1.5}
button{padding:8px 12px;margin:0 8px 12px 0;cursor:pointer}
pre{white-space:pre-wrap;overflow-wrap:anywhere;max-height:320px;overflow:auto}
table{border-collapse:collapse}th,td{text-align:left;padding:8px 16px;border-bottom:1px solid #d8dee6}#alignment{margin-bottom:24px}summary{cursor:pointer}button:focus-visible,summary:focus-visible{outline:3px solid #4477aa}
@media(max-width:800px){.comparison{grid-template-columns:1fr}body{margin:12px}.metrics{min-height:0}}
</style></head>
<body><h1>Src + ligand 4521</h1>
<p>Predicted complexes · Protein: blue cartoon · Ligand: orange sticks</p>
<p>Single-sequence demonstration. Experimental reference IC50: 8.7 nM. Predictions are not experimental measurements.</p>
<section class="panel" id="alignment" hidden>
<h2>Kinase-domain alignment</h2>
<p>OpenFold3: purple · Boltz-2: teal. Both ligands are shown as sticks in their model color.</p>
<p id="alignment-summary"></p>
<button type="button" id="kinase">Kinase domains</button>
<button type="button" id="pockets">Ligands and nearby residues</button>
<button type="button" id="full">Full proteins</button>
<div class="viewer" id="overlay" aria-label="Aligned structure overlay"></div>
<div id="confidence-table"></div>
<details><summary>Alignment method and measurements</summary><pre id="comparison-data"></pre></details>
</section>
<main class="comparison" id="comparison"></main>
<script>
const models=MODELDATA;
const comparison=COMPARISONDATA;
const viewers=[];
let sharedView=null;
const kinaseResidues=Array.from({length:254},(_,i)=>i+270);
if(comparison && typeof $3Dmol!=='undefined'){
    document.getElementById('alignment').hidden=false;
    const c=comparison;
    document.getElementById('alignment-summary').textContent=
        `Aligned on 254 Cα atoms in Src residues 270–523. Kinase RMSD: ${c.kinase_ca_rmsd_angstrom.toFixed(1)} Å. Ligand center separation: ${c.ligand_centroid_distance_after_kinase_fit_angstrom.toFixed(1)} Å. Shared contact residues within 5 Å: ${c.shared_contact_residues.join(', ')}. These predictions disagree substantially even after alignment.`;
    document.getElementById('comparison-data').textContent=JSON.stringify(c,null,2);
    const confidence=c.confidence_from_cif_ca_b_column_0_100;
    document.getElementById('confidence-table').innerHTML=
        '<table><caption>Mean Cα pLDDT (0–100)</caption><thead><tr><th>Region</th><th>OpenFold3</th><th>Boltz-2</th></tr></thead><tbody>'+
        [['Kinase domain','kinase_ca_mean'],['Same residues: union of both ligand pockets','union_pocket_ca_mean']].map(([label,key])=>
            `<tr><td>${label}</td><td>${confidence.openfold3[key].toFixed(1)}</td><td>${confidence.boltz2[key].toFixed(1)}</td></tr>`).join('')+
        '</tbody></table><p>Confidence is model-reported, not experimental validation. Pocket residues have a protein heavy atom within 5 Å of a ligand heavy atom in either prediction.</p>';
    const overlay=$3Dmol.createViewer(document.getElementById('overlay'),{backgroundColor:'white'});
    for(const name of ['openfold3','boltz2']) overlay.addModel(models[name].cif,'cif');
    function showOverlay(mode){
        overlay.setStyle({},{});
        const pocket=[...new Set([...c.contact_residues.openfold3,...c.contact_residues.boltz2])];
        const residues=mode==='pocket'?pocket:kinaseResidues;
        for(const [model,color] of [[0,'#9955bb'],[1,'#008b8b']]){
            const selection=mode==='full'?{model,chain:'A'}:{model,chain:'A',resi:residues};
            overlay.setStyle(selection,{cartoon:{color}});
            if(mode==='pocket') overlay.addStyle(selection,{stick:{color,radius:0.12}});
            overlay.setStyle({model,chain:'B'},{stick:{color,radius:0.25}});
        }
        overlay.zoomTo(mode==='full'?{}:{or:[{chain:'A',resi:residues},{chain:'B'}]});
        overlay.render();
    }
    showOverlay('kinase');
    sharedView=overlay.getView();
    viewers.push(overlay);
    document.getElementById('kinase').onclick=()=>showOverlay('kinase');
    document.getElementById('pockets').onclick=()=>showOverlay('pocket');
    document.getElementById('full').onclick=()=>showOverlay('full');
}
const format=(value,digits=2)=>Number.isFinite(value)?value.toFixed(digits):'Unavailable';
for(const [name,label] of [['openfold3','OpenFold3'],['boltz2','Boltz-2']]){
    const panel=document.createElement('section');
    panel.className='panel';
    panel.innerHTML=`<h2>${label}</h2><p class="metrics"></p>
        <div class="controls"><button type="button" class="whole">Whole complex</button>
        <button type="button" class="ligand">Focus ligand</button></div>
        <div class="viewer" aria-label="${label} predicted complex"></div>
        <details><summary>Model scores</summary><pre></pre></details>`;
    document.getElementById('comparison').append(panel);
    const model=models[name];
    const metrics=panel.querySelector('.metrics');
    if(!model){
        metrics.textContent='No prediction available. Run this model to compare its results.';
        panel.querySelector('.controls').remove();
        panel.querySelector('.viewer').remove();
        panel.querySelector('details').remove();
        continue;
    }
    const a=model.affinity;
    metrics.textContent=a
        ? `Predicted IC50 equivalent: ${format(a.derived_IC50_nM,1)} nM · pIC50: ${format(a.derived_pIC50)} · Binder probability: ${format(a.binder_probability)}`
        : name==='openfold3'
            ? `Complex pLDDT: ${format(model.scores[0]?.complex_plddt_score,1)} / 100. No affinity score.`
            : 'Affinity summary unavailable. See model scores below.';
    panel.querySelector('pre').textContent=JSON.stringify(model.scores,null,2);
    const canvas=panel.querySelector('.viewer');
    if(typeof $3Dmol==='undefined'){
        canvas.textContent='Unable to initialize the bundled 3D library. Regenerate this HTML and reload the preview.';
        panel.querySelectorAll('button').forEach(button=>button.disabled=true);
        continue;
    }
    const viewer=$3Dmol.createViewer(canvas,{backgroundColor:'white'});
    viewer.addModel(model.cif,'cif');
    viewer.setStyle({chain:'A'},{cartoon:{color:'#4477aa'}});
    viewer.setStyle({chain:'B'},{stick:{colorscheme:'orangeCarbon',radius:0.22}});
    viewer.zoomTo();
    if(sharedView) viewer.setView(sharedView);
    viewer.render();
    viewers.push(viewer);
    panel.querySelector('.whole').onclick=()=>{viewer.zoomTo();viewer.render();};
    panel.querySelector('.ligand').onclick=()=>{viewer.zoomTo({chain:'B'});viewer.render();};
}
window.addEventListener('resize',()=>viewers.forEach(viewer=>{viewer.resize();viewer.render();}));
</script></body></html>'''.replace('MODELDATA',data).replace('COMPARISONDATA',json.dumps(comparison)).replace('LIBRARYDATA',library)
(root / 'results' / 'viewer.html').write_text(html)
print('Created results/viewer.html (self-contained; no external script required)')
