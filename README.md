# Src–ligand NVIDIA NIM example

Predict a human Src–ligand complex with OpenFold3 1.5.0 and Boltz-2 1.9.0,
then compare the structures in a browser. Boltz-2 also returns affinity predictions.

[Internal reference: dataset, model inputs, and outputs](docs/overview.html).

## Web app

Start the complete affinity workspace with the existing `.env` and registry login:

```bash
docker compose up --build
```

Open **http://localhost:8080**. Subsequent starts only need `docker compose up`.
The app starts immediately while Boltz-2 warms up and imports the five-ligand
BindingDB example in the background. You can also paste a protein FASTA, add
SMILES, or import a ligand CSV. No manual dataset or Python steps are needed.

- Select a target and multiple ligands; optionally add substitution variants
  such as `T341I` or `T341I,Y530F` (positions are 1-based in the input sequence).
- Run one prediction per ligand–protein pair, then view ranked affinities,
  original/variant comparisons, and interactive complexes with mutation highlights.
- Export CSV/CIF/JSON, reopen saved runs, duplicate a screen, or retry failed jobs.
- `Cancel remaining` cancels queued pairs; the active GPU request finishes.

The app uses a small Python web service and separate prediction worker, plain JavaScript, SQLite, and bundled
3Dmol.js. Its state and results persist in the `app-data` Docker volume; model
caches also persist. `docker compose down` retains these volumes. The app listens
on localhost port 8080 by default. For remote access, forward port 8080 through SSH
or your editor. This is a local workspace without user accounts.

Default GPU assignment: Boltz-2 uses GPU 0 (`BOLTZ_GPU_ID`); the optional OpenFold3
profile uses GPU 1 (`OPENFOLD_GPU_ID`). The web app needs only Boltz-2 and one
supported GPU. `APP_PORT` and `APP_BIND_ADDRESS` configure the published app port.
Run `docker compose --profile openfold up` only when OpenFold3 is also wanted.
Do not share the app's worker endpoint with external benchmark jobs when testing
queue throughput.

Restarting just the web app leaves predictions running in the separate worker.
After an interrupted connection or worker restart during inference, the
queue pauses because the NIM might still be working. Restart Boltz-2, wait for it
to become ready, use **resume the queue**, then **Retry failed** for the interrupted
pair. Completed predictions are retained and previous attempt files are preserved.
See [webapp/README.md](webapp/README.md) for architecture, limits, and testing.

## Long-running BindingDB benchmarks

Open **Benchmarks** in the app to sample thousands of labeled pairs and run them
in the background. The default is 2,000 distinct pairs with exact IC50 labels,
50–1,000 residues, and a fixed seed. Ki and Kd can be selected separately; their
plots are proxy comparisons against predicted IC50, not direct endpoint matches.

Each selected pair retains its source row, label, sequence, SMILES, and provenance.
Results, raw responses, and CIFs are saved after every prediction. Repeated assay
rows do not add duplicate predictions for a sequence–molecule–endpoint pair;
one exact measurement is sampled per pair. The source is the existing September
2026 **BindingDB Articles** archive, not all BindingDB.

The dashboard offers pause/resume, failed-job retry, live scatter plots,
Pearson/Spearman correlations, log-scale errors, and within-target correlations
when enough pairs exist. Export all results as CSV, the reproducible sample as
JSONL, and a standalone SVG plot. Metrics use all successful pairs; the interactive
plot displays at most 5,000 points for responsiveness.

```bash
# Web updates/restarts do not stop the independent inference worker.
docker compose restart app
# Watch the long-running worker.
docker compose logs -f worker
```

The queue and sample live in SQLite and files on the named `app-data` volume.
If sampling is interrupted by a web restart, it is restarted reproducibly from
its saved configuration. Committed predictions are never discarded. A worker
restart marks any in-flight request interrupted and pauses the queue until NIM
is restarted and the operator resumes it; this avoids overlapping an uncertain
old request. Interactive screens take priority between benchmark predictions.

See [benchmark methodology and recovery](webapp/BENCHMARKS.md).

## Script quickstart

Run commands from the repository root. You need:

- A Linux GPU host with Docker Compose v2 and the NVIDIA Container Toolkit
  configured to provide the `nvidia` Docker runtime.
- Two supported NVIDIA GPUs: this Compose configuration assigns GPU 1 to
  OpenFold3 and GPU 0 to Boltz-2. Check the model-specific
  [OpenFold3 hardware requirements](https://docs.nvidia.com/nim/bionemo/openfold3/latest/support-matrix.html)
  and [Boltz-2 support matrix](https://docs.nvidia.com/nim/bionemo/boltz2/latest/support-matrix.html).
- An NGC personal API key with NGC Catalog access and access to both images.
- Python 3.10+ with NumPy for structural alignment, plus internet access
  for dataset and image/model downloads. The generated viewer works offline.

1. Configure credentials. Copy the template only if `.env` does not already exist:

   ```bash
   test -f .env || cp .env.example .env
   chmod 600 .env
   ```

   Edit `.env` to set `NGC_API_KEY`. Authenticate Docker separately:

   ```bash
   docker login nvcr.io --username '$oauthtoken'
   ```

   Paste the same NGC API key at the password prompt. NVIDIA documents this
   authentication in its [getting started guide](https://docs.nvidia.com/nim/bionemo/boltz2/latest/getting-started.html).

2. Start Boltz-2 and run the example (OpenFold3 is disabled by default):

   ```bash
   python3 -m pip install -r requirements.txt
   docker compose up -d
   python3 scripts/run_when_ready.py
   ```

   The runner downloads the pinned September 2026 BindingDB Articles archive
   into `data/raw/` (reusing it on later runs), extracts five records, and
   prepares the model inputs. It then polls readiness every 15 seconds for up to one hour per service,
   runs Boltz-2, and builds `results/viewer.html`.
   First startup downloads model weights; inference can also take time.
   Each prediction has a two-hour HTTP timeout.

   To enable OpenFold3 as well, run `docker compose --profile openfold up -d`
   and `python3 scripts/run_when_ready.py --with-openfold`. Existing OpenFold
   results remain visible in the viewer until removed.

3. Open `results/viewer.html` in your browser. Compare OpenFold3 and Boltz-2 side by side, rotate either complex,
   or choose **Focus ligand** in its panel. On small screens the panels stack. The structures and rendering library are embedded in the HTML.
   In Cursor, use Microsoft Live Preview to open the rendered view.

## Files and outputs

| Path | Purpose |
| --- | --- |
| `compose.yaml` | Pinned NIM images, GPU assignments, local ports, persistent caches |
| `.env.example` | Credential template; `.env` is ignored |
| `data/BindingDB_top5.tsv` | Generated five-record source sample, ignored by Git; the first row drives the example |
| `data/raw/` | Downloaded source archives, ignored by Git |
| `inputs/` | Generated requests, sequence, SMILES, and metadata, ignored by Git |
| `assets/vendor/` | Bundled 3Dmol.js library and licenses for offline viewing |
| `scripts/` | Input preparation, prediction, startup waiting, and viewer generation |
| `results/` | Generated predictions and viewer, ignored by Git |
| [docs/example.md](docs/example.md) | Scientific context and affinity interpretation |

Startup regenerates the sample TSV and `inputs/` from the cached archive.
To prepare data separately, run `python3 scripts/download_bindingdb.py` followed
by `python3 scripts/prepare_example.py`. Predictions write `response.json`, `complex_0.cif`, and `scores.json`
under `results/<model>/`. Boltz-2 also writes `affinity_summary.json`.
Repeated runs overwrite these filenames. The example uses a query-only alignment
and an unverified assay construct; see the interpretation notes before using scores.

## Operate and troubleshoot

```bash
docker compose ps
docker compose logs --tail=100 -f openfold3 boltz2
curl --fail http://localhost:8000/v1/health/ready
curl --fail http://localhost:8001/v1/health/ready
```

A running container may still be loading its model. If startup fails, check logs,
NGC access, GPU visibility (`nvidia-smi`), and available disk/GPU memory.
An HTTP prediction failure saves the server response to `results/<model>/error.txt`.

Once services are ready, individual steps can be rerun:

```bash
python3 scripts/predict.py openfold3
python3 scripts/predict.py boltz2
python3 scripts/make_viewer.py
```

The viewer includes whichever model outputs exist, including outputs from earlier
runs. The runner exits unsuccessfully if any requested prediction fails.

Stop services with `docker compose down`. Named volumes retain model caches and
Boltz-2 server artifacts across restarts. Both APIs bind only to localhost:
OpenFold3 on port 8000 and Boltz-2 on port 8001.

## Structural comparison

When both predictions exist, the viewer generator aligns OpenFold3 onto Boltz-2
using the 254 matching Cα atoms in the [UniProt P12931 kinase domain](https://www.uniprot.org/uniprotkb/P12931/entry)
(residues 270–523). It adds an overlay with domain, pocket, and full-protein views,
and compares confidence over the same pocket residues in both models.
Side-by-side structures initially use the same camera and aligned coordinates;
individual controls then move independently.

Run `python3 scripts/compare_structures.py` to recompute measurements separately.
`results/comparison/summary.json` records the unweighted rigid fit, RMSDs,
ligand center separation, contacts within 5 Å, and mean Cα confidence.
`results/comparison/openfold3_aligned.cif` contains transformed coordinates;
the original predictions are preserved. Ligand atom RMSD is not reported because
atom correspondence and symmetry have not been established. Neither model is
an experimental reference, and a poor domain fit limits precise pose comparison.

## Potency ranking test

For a separate random throughput workload, run `python3 scripts/sample_load_test.py`.
It samples 100 rows without replacement from the cached Articles archive using
seed 20260914, keeping single-chain proteins of 50–1,000 canonical residues with
nonempty ligand SMILES. `inputs/load_test_20260914_100/` contains the source rows,
100 Boltz requests, and a manifest with source row numbers, archive hash, and filters.
Each request uses its row's protein and ligand with the example's inference settings.
SMILES are not chemically validated. This prepares inputs without running inference.
Use `--count`, `--seed`, and `--output` to prepare other samples.

With Boltz-2 ready, run `python3 scripts/rank_potency.py` to predict all five
sample ligands with identical settings. Experimental IC50 labels are used only
for evaluation. Lower IC50 (and lower raw affinity score) means stronger potency.

Each run saves requests, raw responses, `ranking.csv`, `summary.json`, and
`report.md` in a timestamped `results/potency_ranking/` directory. The report
includes strongest-first ranks, Spearman rank correlation, and the number of
concordant, discordant, and tied ligand pairs. This is a small exploratory test;
assay equivalence, training overlap, and run-to-run stability are not verified.
