"""Extract the first BindingDB record and prepare single-sequence NIM requests."""
import csv
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
with (root / 'data' / 'BindingDB_top5.tsv').open(encoding='utf-8', newline='') as f:
    row = next(csv.DictReader(f, delimiter='\t'))
seq = row['BindingDB Target Chain Sequence 1']
smiles = row['Ligand SMILES']
msa = {'main': {'a3m': {'alignment': f'>query\n{seq}\n', 'format': 'a3m'}}}
inputs = root / 'inputs'
inputs.mkdir(exist_ok=True)
(inputs / 'target.fasta').write_text(f'>P12931 Src BindingDB target sequence\n{seq}\n')
(inputs / 'ligand.smi').write_text(f'{smiles} 4521\n')
metadata = {k: row[k] for k in ['BindingDB Reactant_set_id', 'BindingDB MonomerID', 'Target Name', 'IC50 (nM)', 'Article DOI', 'PMID']}
metadata['alignment'] = 'query-only; no homology search'
metadata['target_construct'] = 'Full sequence from BindingDB; assay construct not verified'
(inputs / 'example.json').write_text(json.dumps(metadata, indent=2))
of = {'inputs': [{'input_id': 'src_4521', 'molecules': [
    {'type': 'protein', 'id': 'A', 'sequence': seq, 'msa': msa},
    {'type': 'ligand', 'id': 'B', 'smiles': smiles}], 'diffusion_samples': 1, 'output_format': 'cif'}]}
boltz = {'polymers': [{'id': 'A', 'molecule_type': 'protein', 'sequence': seq, 'msa': msa}],
         'ligands': [{'id': 'B', 'smiles': smiles, 'predict_affinity': True}],
         'recycling_steps': 3, 'sampling_steps': 200, 'diffusion_samples': 1,
         'output_format': 'mmcif', 'sampling_steps_affinity': 200, 'diffusion_samples_affinity': 5}
for name, data in [('openfold3', of), ('boltz2', boltz)]:
    (inputs / f'{name}.json').write_text(json.dumps(data, indent=2))
print(f'Prepared Src ({len(seq)} aa) + ligand 4521')
