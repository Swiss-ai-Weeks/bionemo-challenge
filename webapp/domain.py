"""Validated sequences, mutations, molecules, and the versioned Boltz adapter."""
import math
import re
from rdkit import Chem, RDLogger
from rdkit.Chem import Draw, Descriptors

RDLogger.DisableLog('rdApp.error')
AA = set('ACDEFGHIKLMNPQRSTVWY')
MODEL_VERSION = '1.9.0'
CONVERSION = 'boltz-log10-uM-v1'
DEFAULT_SETTINGS = dict(recycling_steps=3, sampling_steps=200, diffusion_samples=1,
                        sampling_steps_affinity=200, diffusion_samples_affinity=5)


def sequence(text):
    lines = text.strip().splitlines()
    if sum(line.startswith('>') for line in lines) > 1:
        raise ValueError('Provide one protein sequence, not a multi-record FASTA.')
    seq = ''.join(line for line in lines if not line.startswith('>'))
    seq = re.sub(r'\s+', '', seq).upper()
    if not 4 <= len(seq) <= 2048:
        raise ValueError('Protein length must be 4–2048 amino acids for this workspace.')
    invalid = sorted(set(seq) - AA)
    if invalid:
        raise ValueError('Invalid amino-acid codes: ' + ', '.join(invalid))
    return seq


def mutate(parent, changes):
    seq = list(parent)
    seen, normalized = set(), []
    parts = re.split(r'[,;\s]+', changes.strip().upper())
    for part in filter(None, parts):
        match = re.fullmatch(r'([A-Z])(\d+)([A-Z])', part)
        if not match:
            raise ValueError(f'Use substitutions such as T341I; could not read {part}.')
        old, position, new = match.groups()
        position = int(position)
        if old not in AA or new not in AA:
            raise ValueError('Use standard amino-acid codes for substitutions.')
        if position < 1 or position > len(parent):
            raise ValueError(f'Position {position} is outside 1–{len(parent)}.')
        if position in seen:
            raise ValueError(f'Position {position} is specified more than once.')
        if parent[position - 1] != old:
            raise ValueError(f'Position {position} is {parent[position - 1]}, not {old}.')
        if old == new:
            raise ValueError(f'{part} does not change the sequence.')
        seen.add(position)
        seq[position - 1] = new
        normalized.append(dict(position=position, original=old, replacement=new,
                               label=f'{old}{position}{new}',
                               before=parent[max(0, position-7):position+6],
                               after=''.join(seq[max(0, position-7):position+6])))
    if not normalized:
        raise ValueError('Enter at least one substitution, for example T341I.')
    if len(normalized) > 20:
        raise ValueError('Use at most 20 substitutions per variant.')
    result = ''.join(seq)
    for change in normalized:
        p = change['position']
        change['after'] = result[max(0, p-7):p+6]
    return result, sorted(normalized, key=lambda x: x['position'])


def molecule(smiles):
    if not smiles or len(smiles) > 2000:
        raise ValueError('Provide a SMILES string of at most 2000 characters.')
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('The SMILES could not be parsed. Check atom and ring notation.')
    if not 1 <= mol.GetNumHeavyAtoms() <= 128:
        raise ValueError('Ligands must contain 1–128 heavy atoms.')
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError('Provide one connected molecule; remove salts or separate fragments.')
    canonical = Chem.MolToSmiles(mol)
    drawer = Draw.MolDraw2DSVG(180, 110)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return canonical, drawer.GetDrawingText(), round(Descriptors.MolWt(mol), 1)


def settings(data):
    bounds = {'recycling_steps': (1,10), 'sampling_steps': (10,1000),
              'diffusion_samples': (1,5), 'sampling_steps_affinity': (10,1000),
              'diffusion_samples_affinity': (1,10)}
    if set(data) - set(bounds):
        raise ValueError('Unknown prediction setting.')
    result = DEFAULT_SETTINGS | data
    for key, value in result.items():
        low, high = bounds[key]
        if type(value) is not int or not low <= value <= high:
            raise ValueError(f'{key} must be an integer between {low} and {high}.')
    return result


def payload(snapshot):
    seq = snapshot['sequence']
    return dict(polymers=[dict(id='A', molecule_type='protein', sequence=seq,
        msa={'main': {'a3m': {'alignment': f'>query\n{seq}\n', 'format': 'a3m'}}})],
        ligands=[dict(id='B', smiles=snapshot['smiles'], predict_affinity=True)],
        output_format='mmcif', **snapshot['settings'])


def scores(response, model_version=MODEL_VERSION):
    if model_version != MODEL_VERSION:
        raise ValueError('Unsupported model version: affinity conversion must be reviewed.')
    try:
        raw = response['affinities']['B']['affinity_pred_value'][0]
        if isinstance(raw, bool) or not isinstance(raw, (float,int)) or not math.isfinite(raw):
            raise ValueError('Model returned a non-finite affinity score.')
        nm = 1000 * 10 ** raw
        if not math.isfinite(nm) or nm <= 0:
            raise ValueError('Converted affinity is outside the supported numeric range.')
    except (KeyError, IndexError, TypeError, OverflowError) as e:
        raise ValueError('Model response has no valid affinity score for ligand B.') from e
    def optional(value):
        return value if isinstance(value,(int,float)) and math.isfinite(value) else None
    probability = response['affinities']['B'].get('affinity_probability_binary', [None])
    iptm = response.get('ligand_iptm_scores') or [None]
    prob = optional(probability[0]) if probability else None
    if prob is not None and not 0 <= prob <= 1:
        prob = None
    return dict(raw=raw, pic50=6-raw, ic50_nm=nm, binder_probability=prob,
                ligand_iptm=optional(iptm[0]), conversion=CONVERSION)
