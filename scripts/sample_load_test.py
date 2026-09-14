"""Prepare reproducible random BindingDB rows and Boltz requests; no inference."""
import argparse
from collections import Counter
import csv
import hashlib
import io
import json
from pathlib import Path
import random
import zipfile

from download_bindingdb import ARCHIVE, MEMBER

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--count', type=int, default=100)
    parser.add_argument('--seed', type=int, default=20260914)
    parser.add_argument('--min-residues', type=int, default=50)
    parser.add_argument('--max-residues', type=int, default=1000)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    if args.count < 1 or not 1 <= args.min_residues <= args.max_residues:
        parser.error('count must be positive and residue bounds must be ordered and positive')
    archive = ROOT / 'data/raw' / ARCHIVE
    out = args.output or ROOT / 'inputs' / f'load_test_{args.seed}_{args.count}'
    if out.exists():
        parser.error(f'Output already exists: {out}; choose another --output')
    rng = random.Random(args.seed)
    selected = []
    counts = Counter()
    with zipfile.ZipFile(archive) as source, source.open(MEMBER) as member:
        reader = csv.DictReader(io.TextIOWrapper(member, encoding='utf-8'), delimiter='\t')
        fields = reader.fieldnames
        for row_number, row in enumerate(reader, 1):
            counts['total'] += 1
            seq = (row.get('BindingDB Target Chain Sequence 1') or '').strip()
            smiles = (row.get('Ligand SMILES') or '').strip()
            chains = (row.get('Number of Protein Chains in Target (>1 implies a multichain complex)') or '').strip()
            if chains != '1':
                counts['excluded_chain_count'] += 1
                continue
            if not args.min_residues <= len(seq) <= args.max_residues:
                counts['excluded_length'] += 1
                continue
            if set(seq) - set('ACDEFGHIKLMNPQRSTVWY'):
                counts['excluded_sequence_alphabet'] += 1
                continue
            if not smiles:
                counts['excluded_empty_smiles'] += 1
                continue
            counts['eligible'] += 1
            item = (row_number, row, seq, smiles)
            if len(selected) < args.count:
                selected.append(item)
            else:
                index = rng.randrange(counts['eligible'])
                if index < args.count:
                    selected[index] = item
    if len(selected) < args.count:
        raise ValueError(f'Only {len(selected)} eligible rows; requested {args.count}')
    rng.shuffle(selected)
    out.mkdir(parents=True)
    jobs = []
    with (out / 'rows.tsv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        for i, (row_number, row, seq, smiles) in enumerate(selected, 1):
            writer.writerow(row)
            request = {
                'polymers': [{'id': 'A', 'molecule_type': 'protein', 'sequence': seq,
                              'msa': {'main': {'a3m': {'alignment': f'>query\n{seq}\n', 'format': 'a3m'}}}}],
                'ligands': [{'id': 'B', 'smiles': smiles, 'predict_affinity': True}],
                'recycling_steps': 3, 'sampling_steps': 200, 'diffusion_samples': 1,
                'output_format': 'mmcif', 'sampling_steps_affinity': 200,
                'diffusion_samples_affinity': 5,
            }
            name = f'request_{i:04d}.json'
            (out / name).write_text(json.dumps(request, indent=2) + '\n')
            jobs.append({'request': name, 'source_data_row': row_number,
                         'reactant_set_id': row['BindingDB Reactant_set_id'],
                         'ligand_id': row['BindingDB MonomerID'],
                         'target': row['Target Name'], 'protein_residues': len(seq)})
    manifest = {
        'archive': ARCHIVE, 'archive_sha256': hashlib.sha256(archive.read_bytes()).hexdigest(),
        'member': MEMBER, 'seed': args.seed, 'count': args.count,
        'method': 'Uniform reservoir sampling without replacement over eligible source rows, then shuffled.',
        'filters': {'protein_chains': 1, 'min_residues': args.min_residues,
                    'max_residues': args.max_residues, 'canonical_amino_acids_only': True,
                    'nonempty_smiles': True},
        'notes': 'Row sampling retains repeated assay measurements. SMILES are not chemically validated. '
                 'Query-only MSA; source protein sequences are not verified assay constructs. '
                 'This is a throughput workload, not a potency validation dataset.',
        'counts': dict(counts), 'jobs': jobs,
    }
    (out / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(json.dumps({'output': str(out), 'counts': dict(counts), 'selected': len(jobs)}, indent=2))


if __name__ == '__main__':
    main()
