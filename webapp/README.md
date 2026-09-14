# Affinity lab implementation

The Compose `app` service serves a plain HTML/CSS/JavaScript interface and a
FastAPI JSON API. A separate worker container submits one Boltz-2 NIM request at
a time, so web app restarts do not interrupt inference. SQLite holds the library, immutable input snapshots, run queue, and scores;
CIF and raw request/response files live beside it on the persistent `app-data`
volume. The bundled 3Dmol.js library loads only when opening a complex.

## Local use

With NVIDIA registry access and `NGC_API_KEY` already configured, run
`docker compose up --build` from the repository root and open localhost:8080.
Only one GPU is required. On later starts `docker compose up` is sufficient.
A background download imports the pinned BindingDB five-ligand example with MD5
verification. If it fails, use the UI retry button or enter custom inputs.

Protein sequences accept standard amino acids, one chain, and 4–2048 residues.
Variants support up to 20 simultaneous substitutions, validated against the
immutable parent sequence. No insertions/deletions or PDB numbering translation.
Ligands accept one connected SMILES with up to 128 heavy atoms. RDKit validates
molecules, canonicalizes SMILES for deduplication, and renders SVG thumbnails.

CSV columns are `name,smiles`, optionally `endpoint,value,unit,source`. Reference
values must be in nM; qualifiers are retained. Partial import reports errors by
row. A reference belongs to its original target and is not copied to variants.

Runs permit up to 50 ligands, 6 protein sequences, and 100 total pairs. Each pair
gets one ligand and one protein chain. Scores rank within a protein variant; exact
ties share an average rank. Missing/failed scores stay unavailable. The ranking table displays ground-truth measurements alongside predictions,
preserving their endpoint, unit, source and inequality qualifier. The comparison
table includes the original protein ground truth and shows pIC50 and its difference from the original sequence, and CSV includes
both absolute predictions and differences. Each run uses identical settings for
all pairs and stores versioned provenance. Query-only alignments are explicitly
labeled; small score differences and mutation effects are exploratory.

NIM 1.9.0 calls `(6 - raw) * 1.364` `affinity_pic50`. We preserve this raw response
but derive displayed values from upstream Boltz's log10(micromolar) convention:
`pIC50 = 6 - raw`; `IC50 equivalent in nM = 1000 * 10 ** raw`. The conversion is
versioned in `domain.py` and checked in tests. Binder probability and interface
confidence are separate quantities.

## Recovery and operations

- Exactly one inference worker is supported, enforced with a shared file lock.
  Do not scale worker replicas; keep the web service at one process. SQLite connections use short transactions.
- Queued jobs persist through restarts. On an inference worker restart, active jobs become interrupted; the queue
  pauses until the operator restarts NIM and explicitly resumes it. Readiness
  alone cannot prove an old inference stopped.
- Connection timeouts and NIM 5xx responses also pause the queue. Definitive 4xx
  failures retain other completed results and allow subsequent pairs to run.
- Retries keep previous attempt files and diagnostics. `Cancel remaining` cancels
  queued jobs, allowing the active request to finish. Closing a tab has no effect.
- Prediction timeout is two hours; the readiness monitor checks every five seconds.
- Run history displays the latest 100 screens. There is no automatic deletion.
- Raw files are downloadable through allowlisted artifact kinds, not filesystem
  paths. CSV labels are escaped against spreadsheet formulas. Browser writes are
  same-origin; no credentials are exposed to the app UI.
- The default localhost binding is intentional. Add authentication and access
  control before making this multi-user or internet-facing.

## Development and tests

Use Python 3.12. Install `webapp/requirements.txt`, plus pytest and httpx for tests:

```bash
python -m pip install -r webapp/requirements.txt pytest==9.1.1 httpx==0.28.1
python -m pytest webapp/tests -q
APP_DATA_DIR=/tmp/affinity-dev BOLTZ_URL=http://localhost:8001 \
  uvicorn webapp.main:app --host 127.0.0.1 --port 8080
```

Tests use temporary databases and mocked NIM responses; no GPU is needed. They
cover validation, input snapshots, reference isolation, affinity conversion, ties,
cancellation, failed-job isolation, retries, restart recovery, and CSV imports.
The static files have no frontend build step. Rebuild the app image after changes:
`docker compose build app worker`, then restart only idle services. Web-only
restarts are safe during predictions: `docker compose up -d --no-deps app`.

The existing standalone analysis/download scripts remain usable. The new app
reuses the archive downloader but maintains its own database and artifacts.

See [BENCHMARKS.md](BENCHMARKS.md) for long-running random sampling and analysis.
