# Src–ligand NVIDIA NIM example

Predict a human Src–ligand complex with OpenFold3 1.5.0 and Boltz-2 1.9.0,
then compare the structures in a browser. Boltz-2 also returns affinity predictions.

[Internal reference: dataset, model inputs, and outputs](docs/overview.html).

## Quickstart

Run commands from the repository root. You need:

- A Linux GPU host with Docker Compose v2 and the NVIDIA Container Toolkit
  configured to provide the `nvidia` Docker runtime.
- Two supported NVIDIA GPUs: this Compose configuration assigns GPU 0 to
  OpenFold3 and GPU 1 to Boltz-2. Check the model-specific
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

2. Start the services and run the example:

   ```bash
   python3 -m pip install -r requirements.txt
   docker compose up -d
   python3 scripts/run_when_ready.py
   ```

   The runner downloads the pinned September 2026 BindingDB Articles archive
   into `data/raw/` (reusing it on later runs), extracts five records, and
   prepares the model inputs. It then polls readiness every 15 seconds for up to one hour per service,
   runs both predictions concurrently, and builds `results/viewer.html`.
   First startup downloads model weights; inference can also take time.
   Each prediction has a two-hour HTTP timeout.

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
runs. The runner exits unsuccessfully if either prediction fails.

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

With Boltz-2 ready, run `python3 scripts/rank_potency.py` to predict all five
sample ligands with identical settings. Experimental IC50 labels are used only
for evaluation. Lower IC50 (and lower raw affinity score) means stronger potency.

Each run saves requests, raw responses, `ranking.csv`, `summary.json`, and
`report.md` in a timestamped `results/potency_ranking/` directory. The report
includes strongest-first ranks, Spearman rank correlation, and the number of
concordant, discordant, and tied ligand pairs. This is a small exploratory test;
assay equivalence, training overlap, and run-to-run stability are not verified.
