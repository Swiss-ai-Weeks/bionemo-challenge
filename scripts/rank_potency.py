"""Run the five BindingDB example ligands through Boltz-2 and assess potency order."""
import copy
import csv
from datetime import datetime, timezone
import itertools
import json
import math
from pathlib import Path
import statistics
import urllib.error
import urllib.request

ROOT = Path(__file__).resolve().parents[1]
ENDPOINT = 'http://127.0.0.1:8001/biology/mit/boltz2/predict'


def ranks(values):
    """Ascending ranks, averaging tied positions. Lower IC50 = stronger potency."""
    ordered = sorted(values)
    return [statistics.mean(i + 1 for i, other in enumerate(ordered) if other == value)
            for value in values]


def evaluate(records):
    experimental = [r['experimental_ic50_nm'] for r in records]
    predicted = [r['raw_affinity_pred_value'] for r in records]
    er, pr = ranks(experimental), ranks(predicted)
    rho = statistics.correlation(er, pr) if len(set(pr)) > 1 and len(set(er)) > 1 else None
    pairs = []
    for i, j in itertools.combinations(range(len(records)), 2):
        e, p = experimental[i] - experimental[j], predicted[i] - predicted[j]
        status = 'tied' if e == 0 or p == 0 else 'concordant' if e * p > 0 else 'discordant'
        pairs.append({'ligands': [records[i]['ligand_id'], records[j]['ligand_id']], 'status': status})
    for row, e, p in zip(records, er, pr):
        row.update(experimental_rank=e, predicted_rank=p)
    return {
        'spearman_rho': rho,
        'exact_order_preserved': er == pr,
        'concordant_pairs': sum(p['status'] == 'concordant' for p in pairs),
        'discordant_pairs': sum(p['status'] == 'discordant' for p in pairs),
        'tied_pairs': sum(p['status'] == 'tied' for p in pairs),
        'total_pairs': len(pairs),
        'experimental_order_strongest_first': [r['ligand_id'] for r in sorted(records, key=lambda r:r['experimental_rank'])],
        'predicted_order_strongest_first': [r['ligand_id'] for r in sorted(records, key=lambda r:r['predicted_rank'])],
        'records': records, 'pairs': pairs,
    }


def main():
    with (ROOT / 'data/BindingDB_top5.tsv').open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    template = json.loads((ROOT / 'inputs/boltz2.json').read_text())
    if len(rows) != 5 or len({r['BindingDB MonomerID'] for r in rows}) != 5:
        raise ValueError('Expected five distinct example ligands')
    if any(r['BindingDB Target Chain Sequence 1'] != template['polymers'][0]['sequence'] for r in rows):
        raise ValueError('Target sequences differ from the prepared Src request')
    if len({r['Article DOI'] for r in rows}) != 1:
        raise ValueError('Expected measurements from one source paper')
    # Reject censored measurements rather than turning bounds into exact ranks.
    labels = [float(r['IC50 (nM)'].strip()) for r in rows]
    if any(not math.isfinite(v) or v <= 0 for v in labels):
        raise ValueError('Expected positive, finite IC50 values')
    out = ROOT / 'results/potency_ranking' / datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out.mkdir(parents=True)
    records = []
    for i, (row, label) in enumerate(zip(rows, labels), 1):
        ligand = row['BindingDB MonomerID']
        folder = out / ligand
        folder.mkdir()
        request = copy.deepcopy(template)
        request['ligands'][0]['smiles'] = row['Ligand SMILES']
        (folder / 'request.json').write_text(json.dumps(request, indent=2))
        print(f'[{i}/5] Predicting ligand {ligand}', flush=True)
        req = urllib.request.Request(ENDPOINT, data=json.dumps(request).encode(), headers={'Content-Type': 'application/json'})
        try:
            with urllib.request.urlopen(req, timeout=7200) as response:
                result = json.load(response)
        except urllib.error.HTTPError as error:
            (folder / 'error.txt').write_bytes(error.read())
            raise
        (folder / 'response.json').write_text(json.dumps(result, indent=2))
        raw = result['affinities']['B']['affinity_pred_value'][0]
        if not math.isfinite(raw):
            raise ValueError(f'Non-finite affinity for {ligand}')
        records.append({
            'ligand_id': ligand, 'reactant_set_id': row['BindingDB Reactant_set_id'],
            'article_doi': row['Article DOI'], 'experimental_ic50_nm': label,
            'raw_affinity_pred_value': raw, 'derived_predicted_ic50_nm': 1000 * 10 ** raw,
            'binder_probability': result['affinities']['B']['affinity_probability_binary'][0],
        })
        print(f'  Experimental: {label:g} nM; predicted: {1000 * 10 ** raw:.2f} nM', flush=True)
    summary = evaluate(records)
    summary['method'] = 'One fresh Boltz-2 NIM 1.9.0 request per ligand, identical template settings; only ligand SMILES changed. Experimental labels were not submitted. Ascending raw affinity score ranks strongest first; IC50 conversion is monotonic and does not affect ranking. Ties receive average ranks.'
    summary['limitations'] = 'Five same-paper ligands; assay-condition equivalence and training overlap not verified. One request per ligand does not measure run-to-run ranking stability.'
    (out / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    with (out / 'ranking.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(sorted(records, key=lambda r:r['experimental_rank']))
    lines = ['# Boltz-2 potency ranking', '', summary['method'], '',
             f"Spearman rho: {summary['spearman_rho']}; exact order preserved: {summary['exact_order_preserved']}.",
             f"Concordant pairs: {summary['concordant_pairs']}/{summary['total_pairs']}; discordant: {summary['discordant_pairs']}; tied: {summary['tied_pairs']}.", '',
             '| Ligand | Experimental IC50 (nM) | Predicted IC50 (nM) | Experimental rank | Predicted rank |',
             '| --- | ---: | ---: | ---: | ---: |']
    for r in sorted(records, key=lambda r:r['experimental_rank']):
        lines.append(f"| {r['ligand_id']} | {r['experimental_ic50_nm']:g} | {r['derived_predicted_ic50_nm']:.2f} | {r['experimental_rank']:g} | {r['predicted_rank']:g} |")
    lines += ['', summary['limitations'], '']
    (out / 'report.md').write_text('\n'.join(lines))
    print('\n'.join(lines), flush=True)
    print(f'Saved run to {out}', flush=True)


if __name__ == '__main__':
    main()
