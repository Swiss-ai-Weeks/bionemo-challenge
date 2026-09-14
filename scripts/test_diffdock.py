"""Dock ligand 4521 onto the existing Boltz-2 Src structure; save poses and latency."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time
import urllib.error
import urllib.request
from compare_structures import read_atoms

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--url', default='http://127.0.0.1:8003')
    parser.add_argument('--poses', type=int, default=10)
    parser.add_argument('--runs', type=int, default=2)
    args = parser.parse_args()
    if args.poses < 1 or args.runs < 1:
        parser.error('poses and runs must be positive')
    source = ROOT / 'results/boltz2/complex_0.cif'
    atoms = [a for a in read_atoms(source)[2] if a['chain'] == 'A']
    if len({a['residue'] for a in atoms}) != 536:
        raise ValueError('Expected full 536-residue Src protein')
    lines=[]
    for i,a in enumerate(atoms,1):
        name=a['name'] if len(a['name'])==4 else ' '+a['name']
        x,y,z=a['xyz']
        lines.append(f"ATOM  {i:5d} {name:<4} {a['residue_name']:>3} A{a['residue']:4d}    {x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{a['confidence']:6.2f}          {a['element']:>2}  ")
    protein='\n'.join(lines)+'\nTER\nEND\n'
    smiles=json.loads((ROOT/'inputs/boltz2.json').read_text())['ligands'][0]['smiles']
    request={'protein':protein,'ligand':smiles+'\n','ligand_file_type':'txt',
             'num_poses':args.poses,'time_divisions':20,'steps':18,
             'save_trajectory':False,'is_staged':False}
    out=ROOT/'results/diffdock'/datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True)
    (out/'target.pdb').write_text(protein)
    (out/'request.json').write_text(json.dumps(request,indent=2))
    summary={'image':'nvcr.io/nim/mit/diffdock:2.3.0',
             'receptor_source':str(source.relative_to(ROOT)),
             'target':'Src P12931, 536 residues, Boltz-2 predicted receptor; original ligand removed',
             'ligand_id':'4521','num_poses':args.poses,'runs':[],
             'note':'Pose confidence is not affinity. This tests docking on a predicted receptor, not experimental pose accuracy. HTTP times include preprocessing and response transfer, exclude deployment and protein folding.'}
    for i in range(args.runs):
        folder=out/f'run_{i+1}'
        folder.mkdir()
        print(f'Run {i+1}/{args.runs}: requesting {args.poses} poses',flush=True)
        start=time.perf_counter()
        req=urllib.request.Request(args.url+'/molecular-docking/diffdock/generate',data=json.dumps(request).encode(),headers={'Content-Type':'application/json'})
        try:
            with urllib.request.urlopen(req,timeout=1800) as response:
                result=json.load(response)
        except urllib.error.HTTPError as error:
            (folder/'error.txt').write_bytes(error.read())
            raise
        elapsed=time.perf_counter()-start
        (folder/'response.json').write_text(json.dumps(result,indent=2))
        status=result.get('status')
        if status!='success' and status!=['success']:
            raise RuntimeError(f'DiffDock failed: {result.get("details")} (see {folder})')
        poses=result['ligand_positions']; confidence=result['position_confidence']
        if poses and isinstance(poses[0],list):
            poses,confidence=poses[0],confidence[0]
        if len(poses)!=args.poses or len(confidence)!=args.poses:
            raise ValueError('Unexpected number of poses or confidence scores')
        ranked=sorted(zip(confidence,poses),key=lambda p:p[0],reverse=True)
        for rank,(score,sdf) in enumerate(ranked,1):
            if 'M  END' not in sdf:
                raise ValueError('Invalid SDF pose')
            (folder/f'pose_{rank:02d}.sdf').write_text(sdf)
        summary['runs'].append({'run':i+1,'http_seconds':elapsed,'poses_per_second':len(poses)/elapsed,
                                'confidence_descending':[c for c,p in ranked]})
        (out/'summary.json').write_text(json.dumps(summary,indent=2))
        print(json.dumps(summary['runs'][-1]),flush=True)
    print(f'Saved {out}',flush=True)


if __name__=='__main__':
    main()
