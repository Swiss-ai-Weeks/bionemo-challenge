"""Regenerate static documentation assets; optional dependencies: RDKit and NumPy."""
import json
from pathlib import Path
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from compare_structures import read_atoms

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'docs/assets'
OUT.mkdir(exist_ok=True)
smiles = json.loads((ROOT / 'inputs/boltz2.json').read_text())['ligands'][0]['smiles']
molecule = Chem.MolFromSmiles(smiles)
if molecule is None:
    raise ValueError('Invalid ligand SMILES')
rdDepictor.Compute2DCoords(molecule)
drawer = rdMolDraw2D.MolDraw2DSVG(700, 310)
drawer.DrawMolecule(molecule)
drawer.FinishDrawing()
(OUT / 'ligand-4521.svg').write_text(drawer.GetDrawingText())

# Orthographic projection of actual predicted C-alpha coordinates, not a fold cartoon.
atoms = read_atoms(ROOT / 'results/boltz2/complex_0.cif')[2]
ca = [a for a in atoms if a['chain'] == 'A' and a['name'] == 'CA']
xyz = np.array([a['xyz'] for a in ca])
center = xyz.mean(axis=0)
_, _, axes = np.linalg.svd(xyz-center, full_matrices=False)
projected = (xyz-center) @ axes[:2].T
scale = min(510 / np.ptp(projected[:,0]), 260 / np.ptp(projected[:,1]))
def project(points):
    p = (np.array(points)-center) @ axes[:2].T * scale
    return p * [1,-1] + [340,155]
xy = project(xyz)
svg=['<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 700 340" role="img"><title>Boltz-2 predicted Src structure, C-alpha trace</title><rect width="700" height="340" fill="white"/>']
for i in range(1,len(ca)):
    r=ca[i]['residue']
    color = '#9955bb' if 84<=r<=145 else '#2775b6' if 151<=r<=248 else '#008b8b' if 270<=r<=523 else '#b0b8c3'
    svg.append(f'<path d="M{xy[i-1,0]:.2f},{xy[i-1,1]:.2f} L{xy[i,0]:.2f},{xy[i,1]:.2f}" fill="none" stroke="{color}" stroke-width="2"/>')
ligand=np.mean([a['xyz'] for a in atoms if a['chain']=='B' and a['element']!='H'],axis=0)
x,y=project([ligand])[0]
svg.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="6" fill="#df8a22" stroke="white"/><text x="{x+10:.2f}" y="{y:.2f}" font-family="system-ui" font-size="13">Ligand center</text>')
svg.append('<text x="25" y="328" font-family="system-ui" font-size="12" fill="#566273">Predicted coordinates · 2D projection · not an experimental structure</text></svg>')
(OUT / 'src-predicted-trace.svg').write_text(''.join(svg))
print('Created ligand graph and predicted protein trace in docs/assets/')
