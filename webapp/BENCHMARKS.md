# Persistent BindingDB benchmarks

A benchmark is a saved sampling specification plus a persistent set of jobs. It
uses the existing pinned, checksum-verified BindingDB Articles September 2026 ZIP.
This subset is suitable for an exploratory comparison, not evidence of a held-out
benchmark: training overlap and exact experimental constructs are unverified.

## Sampling

- Choose a count (default 2,000, maximum 50,000), fixed seed, endpoints, and residue
  limits (default 50–1,000). Only single-chain standard-amino-acid sequences pass.
- Labels must be exact, finite, positive values in nM. Inequalities and approximate
  values are excluded and counted. IC50, Ki, and Kd are never pooled into one plot.
- Ligands must parse in RDKit, be connected, and have 1–128 heavy atoms. Canonical
  isomeric SMILES defines chemical identity; stereochemistry is retained.
- Uniform reservoir sampling selects distinct `(sequence, canonical SMILES,
  endpoint)` pairs without replacement. Sampling is not uniform over targets.
- Repeated source measurements for a selected pair are sampled uniformly to retain
  one measurement, not averaged or counted as independent predictions. A separate
  seeded random stream chooses the representative label.
- Multiple selected endpoints can produce separate jobs for the same molecule and
  sequence; each is an independent endpoint-specific comparison.
- The sample, source row numbers, IDs, labels, source dates, sequences, and SMILES
  are written to `benchmarks/<id>/sample.jsonl` in the app data volume. The manifest
  stores filtering counts, archive SHA-256, seed, settings, and methodology.
- If there are too few eligible distinct pairs, preparation fails clearly without
  submitting a smaller or duplicated sample. Create a new specification instead.

All jobs use the versioned Boltz-2 NIM 1.9.0 conversion and default query-only
inference settings. Experimental labels never appear in the prediction request.
Ordinary web screens are scheduled ahead of benchmark jobs, between predictions.

## Durability and restart behavior

The web service and inference worker share the named `app-data` volume. SQLite
uses WAL transactions. The worker has an exclusive file lock and the transactional
claim also prevents a second running job. Each request has a stable job ID and
numbered attempt directories with request, response, structure, and diagnostics.

- Browser refresh/close: no effect on sampling or inference.
- `docker compose restart app`: running inference continues in `worker`.
- App restart while sampling: configuration is retained; sampling restarts with
  the same seed. Jobs are inserted in one transaction only after selection is
  complete, so sampling recovery does not partially or doubly enqueue a sample.
- Worker crash/restart: completed jobs remain saved. An in-flight job becomes
  interrupted; the global queue pauses because client loss cannot prove NIM
  cancelled its computation. Restart NIM, wait for readiness, resume the queue,
  and retry interrupted/failed jobs. Pending work is retained.
- Pause benchmark: no new benchmark pairs are claimed; the active request finishes.
- Cancel remaining: pending pairs become cancelled; historical outputs remain.
- `docker compose down`: named volumes remain. Deleting volumes (`down -v`) removes
  the workspace, so do not use it to stop an ongoing benchmark.

## Analysis

Each successful pair supplies experimental `pX = 9 - log10(label_nM)` and predicted
`pIC50 = 6 - raw_score`. IC50 is a matched-endpoint comparison; Ki/Kd versus predicted
IC50 are explicitly proxy comparisons. Scatter plots use equal x/y limits and an
identity line. Missing and failed predictions are excluded from metrics but remain
in the exported records and completion counts.

Pearson and Spearman need at least three points and nonconstant values. Exact ties
receive average ranks. MAE/RMSE are on the log10 molar scale, not raw nM. Pooled
correlation across targets can be driven by between-target differences, so the
report also lists within-target correlations for targets with at least five
completed pairs. No confidence intervals or ranking-reliability claim is inferred
from a single prediction per pair.

The interactive view plots up to 5,000 uniformly selected successful points using
a fixed display seed. Metrics always use all successful pairs for the chosen
endpoint and target filter. The downloadable SVG is rendered by Matplotlib. CSV
contains all jobs (including failures), source identities, ground truth, predicted
scores, confidence values, timestamps, and attempt numbers.

## API

- `POST /api/benchmarks`: name, count, seed, endpoints, min/max residues, settings.
- `GET /api/benchmarks`: summaries; `GET /api/benchmarks/<id>`: progress.
- `POST /api/benchmarks/<id>/pause|resume|cancel|retry`.
- `GET /api/benchmarks/<id>/report?endpoint=IC50&target=...`: statistics and points.
- `GET /api/benchmarks/<id>/records?offset=0&limit=50`: paginated saved pairs.
- `GET /api/benchmarks/<id>/export|manifest|sample|plot.svg`: downloadable artifacts.

API and UI changes do not require restarting the inference worker. After changing
shared worker code, pause benchmark scheduling and allow the active prediction
to finish before deploying the new worker image.
