"""Compare generated Src complexes using a kinase-domain C-alpha rigid fit.

The atom reader handles the single-line atom-site rows emitted by these NIMs;
it deliberately rejects unsupported/malformed rows rather than guessing.
"""
import json
from pathlib import Path
import re
import shlex

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DOMAIN = (270, 523)  # UniProt P12931 protein kinase domain, inclusive.


def read_atoms(path):
    lines = path.read_text().splitlines()
    headers, atoms = [], []
    for line_number, line in enumerate(lines):
        if line.startswith('_atom_site.'):
            headers.append(line.strip().split('.', 1)[1])
        elif headers and line.startswith(('ATOM ', 'HETATM ')):
            values = shlex.split(line)
            if len(values) != len(headers):
                raise ValueError(f'Unsupported atom-site row in {path}:{line_number + 1}')
            row = dict(zip(headers, values))
            if row['pdbx_PDB_model_num'] != '1':
                raise ValueError('Expected one structure model')
            atoms.append(dict(
                line=line_number, values=values, chain=row['label_asym_id'],
                residue=int(row['label_seq_id']) if row['label_seq_id'].isdigit() else None,
                residue_name=row['label_comp_id'], name=row['label_atom_id'],
                element=row['type_symbol'], confidence=float(row['B_iso_or_equiv']),
                xyz=np.array([float(row[f'Cartn_{axis}']) for axis in 'xyz']),
            ))
    if not atoms:
        raise ValueError(f'No atoms in {path}')
    return lines, headers, atoms


def fit(mobile, reference):
    """Least-squares proper rotation for row-vector coordinates (no reflection)."""
    x, y = mobile.mean(axis=0), reference.mean(axis=0)
    u, _, vt = np.linalg.svd((mobile - x).T @ (reference - y))
    correction = np.diag([1, 1, np.linalg.det(u @ vt)])
    rotation = u @ correction @ vt
    return rotation, y - x @ rotation


def rmsd(a, b):
    return float(np.sqrt(np.mean(np.sum((a - b) ** 2, axis=1))))


def compare(results=ROOT / 'results'):
    parsed = {m: read_atoms(results / m / 'complex_0.cif') for m in ['openfold3', 'boltz2']}
    atoms = {m: p[2] for m, p in parsed.items()}
    ca = {m: {a['residue']: a for a in aa if a['chain'] == 'A' and a['name'] == 'CA'}
          for m, aa in atoms.items()}
    if any(set(a) != set(range(1, 537)) for a in ca.values()):
        raise ValueError('Comparison requires the complete 536-residue Src example')
    if any(ca['openfold3'][r]['residue_name'] != ca['boltz2'][r]['residue_name'] for r in ca['openfold3']):
        raise ValueError('Protein residue identities do not match')
    domain = list(range(DOMAIN[0], DOMAIN[1] + 1))
    coords = lambda m, residues: np.array([ca[m][r]['xyz'] for r in residues])
    rotation, translation = fit(coords('openfold3', domain), coords('boltz2', domain))
    aligned = lambda xyz: xyz @ rotation + translation
    out = results / 'comparison'
    out.mkdir(exist_ok=True)
    lines, headers, aa = parsed['openfold3']
    for atom in aa:
        # Preserve all original CIF tokens except Cartesian coordinates.
        tokens = re.findall(r"'(?:[^']*)'|\"(?:[^\"]*)\"|\S+", lines[atom['line']])
        if len(tokens) != len(headers):
            raise ValueError('Unsupported quoted atom row')
        for axis, value in zip('xyz', aligned(atom['xyz'])):
            tokens[headers.index(f'Cartn_{axis}')] = f'{value:.6f}'
        lines[atom['line']] = ' '.join(tokens)
    (out / 'openfold3_aligned.cif').write_text('\n'.join(lines) + '\n')
    ligands = {m: [a for a in aa if a['chain'] == 'B' and a['element'] not in ['H', 'D']]
               for m, aa in atoms.items()}
    if any(not a for a in ligands.values()):
        raise ValueError('Missing ligand heavy atoms')
    ligand_coords = {m: np.array([a['xyz'] for a in aa]) for m, aa in ligands.items()}
    contacts = {}
    confidence = {}
    for m, aa in atoms.items():
        protein = [a for a in aa if a['chain'] == 'A' and a['element'] not in ['H', 'D']]
        distances = np.linalg.norm(np.array([a['xyz'] for a in protein])[:, None, :] - ligand_coords[m][None, :, :], axis=2)
        contacts[m] = sorted({a['residue'] for a, d in zip(protein, distances.min(axis=1)) if d <= 5})
    pocket = sorted(set(contacts['openfold3']) | set(contacts['boltz2']))
    for m in atoms:
        confidence[m] = {
            'kinase_ca_mean': float(np.mean([ca[m][r]['confidence'] for r in domain])),
            'own_pocket_ca_mean': float(np.mean([ca[m][r]['confidence'] for r in contacts[m]])) if contacts[m] else None,
            'union_pocket_ca_mean': float(np.mean([ca[m][r]['confidence'] for r in pocket])) if pocket else None,
        }
    summary = {
        'method': 'Unweighted Kabsch fit of OpenFold3 onto Boltz-2 using corresponding kinase C-alpha atoms; no outlier rejection.',
        'domain_source': 'https://www.uniprot.org/uniprotkb/P12931/entry',
        'kinase_residue_range': list(DOMAIN), 'fit_atom_count': len(domain),
        'kinase_ca_rmsd_angstrom': rmsd(aligned(coords('openfold3', domain)), coords('boltz2', domain)),
        'whole_protein_ca_rmsd_after_kinase_fit_angstrom': rmsd(aligned(coords('openfold3', range(1, 537))), coords('boltz2', range(1, 537))),
        'ligand_centroid_distance_after_kinase_fit_angstrom': float(np.linalg.norm(aligned(ligand_coords['openfold3']).mean(axis=0) - ligand_coords['boltz2'].mean(axis=0))),
        'contact_cutoff_angstrom': 5, 'contact_residues': contacts,
        'shared_contact_residues': sorted(set(contacts['openfold3']) & set(contacts['boltz2'])),
        'confidence_from_cif_ca_b_column_0_100': confidence,
        'confidence_note': 'C-alpha pLDDT from generated CIF B columns; model confidence is not experimental validation.',
        'ligand_note': 'Centroid separation and residue contacts describe site agreement; ligand atom RMSD is not computed because atom correspondence/symmetry has not been established.',
        'pocket_ca_rmsd_after_kinase_fit_angstrom': rmsd(aligned(coords('openfold3', pocket)), coords('boltz2', pocket)) if pocket else None,
        'rotation_row_vectors': rotation.tolist(), 'translation': translation.tolist(),
    }
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    return summary


if __name__ == '__main__':
    print(json.dumps(compare(), indent=2))
