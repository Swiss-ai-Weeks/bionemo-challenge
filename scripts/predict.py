"""Run a prepared request against a local NVIDIA NIM and save its artifacts."""
import argparse
import json
from pathlib import Path
import urllib.request
import urllib.error

root = Path(__file__).resolve().parents[1]
p = argparse.ArgumentParser()
p.add_argument('model', choices=['openfold3', 'boltz2'])
args = p.parse_args()
model = args.model
port, endpoint = (8000, 'openfold/openfold3') if model == 'openfold3' else (8001, 'mit/boltz2')
request = urllib.request.Request(f'http://127.0.0.1:{port}/biology/{endpoint}/predict',
    data=(root / 'inputs' / f'{model}.json').read_bytes(), headers={'Content-Type': 'application/json'})
print(f'Running {model} on Src + ligand 4521...', flush=True)
out = root / 'results' / model
out.mkdir(parents=True, exist_ok=True)
try:
    with urllib.request.urlopen(request, timeout=7200) as response:
        result = json.load(response)
except urllib.error.HTTPError as error:
    detail = error.read().decode('utf-8', errors='replace')
    (out / 'error.txt').write_text(detail)
    raise RuntimeError(f'{model}: HTTP {error.code}: {detail}') from error
(out / 'response.json').write_text(json.dumps(result, indent=2))
if model == 'openfold3':
    structures = result['outputs'][0]['structures_with_scores']
    scores = [{k: v for k, v in s.items() if k != 'structure'} for s in structures]
else:
    structures = result['structures']
    scores = {k: v for k, v in result.items() if k != 'structures'}
if not structures:
    raise RuntimeError('No predicted structures returned')
for i, structure in enumerate(structures):
    (out / f'complex_{i}.cif').write_text(structure['structure'])
(out / 'scores.json').write_text(json.dumps(scores, indent=2))
if model == 'boltz2':
    raw = result['affinities']['B']['affinity_pred_value'][0]
    summary = {
        'raw_log10_ic50_uM': raw,
        'derived_pIC50': 6 - raw,
        'derived_IC50_nM': 1000 * 10 ** raw,
        'experimental_IC50_nM': 8.7,
        'binder_probability': result['affinities']['B']['affinity_probability_binary'][0],
        'note': 'NIM 1.9.0 affinity_pic50 uses (6 - raw) * 1.364; do not interpret that field as dimensionless pIC50. Derived values here use the upstream Boltz log10(IC50 in micromolar) convention.'
    }
    (out / 'affinity_summary.json').write_text(json.dumps(summary, indent=2))
print(json.dumps(scores, indent=2), flush=True)
print(f'Saved artifacts in {out}', flush=True)
