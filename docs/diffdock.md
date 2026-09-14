# Local DiffDock test

DiffDock NIM 2.3.0 docks ligand 4521 onto the existing Boltz-2 predicted Src
receptor (536 residues). Chain B is removed; chain A is exported as PDB.
This checks local inference and candidate pose generation, not experimental accuracy.

```bash
docker compose -f compose.diffdock.yaml up -d
curl --fail http://127.0.0.1:8003/v1/health/ready
python3 scripts/test_diffdock.py --poses 10 --runs 2
python3 scripts/make_diffdock_viewer.py
```

Wait for readiness before running the test. The first startup downloads weights.
The separate Compose project uses GPU 0 and port 8003 by default; override with
`DIFFDOCK_GPU_ID` and `DIFFDOCK_PORT`. It reuses the existing `.env` NGC key.
Model weights persist in its named volume.

Each timestamped `results/diffdock/` folder contains `target.pdb`, the JSON
request, per-run raw responses and confidence-ranked SDF poses, and `summary.json`.
The viewer generator writes a self-contained `viewer.html` for the latest run.
Timings measure the full HTTP call, including preprocessing and response transfer;
they exclude image/model downloads and receptor folding. Pose confidence is not
a binding-affinity score. The test uses SMILES input, 20 time divisions, 18 steps,
and 10 poses by default.

```bash
docker compose -f compose.diffdock.yaml logs --tail=100 diffdock
docker compose -f compose.diffdock.yaml down
```

Sources: [NIM setup](https://docs.nvidia.com/nim/bionemo/diffdock/latest/getting-started.html),
[SMILES and pose sampling](https://docs.nvidia.com/nim/bionemo/diffdock/latest/advanced-usage.html).

## Local result · 14 September 2026

H100 NVL, Src + ligand 4521, 10 poses per request:

| Request | Full HTTP latency | Poses/second |
| --- | ---: | ---: |
| First after startup | 2.95 s | 3.39 |
| Repeat | 1.14 s | 8.74 |

Both requests succeeded. All 20 SDFs parsed as 3D molecules with the input ligand's
chemical identity. These two timings are a smoke test, not a throughput benchmark.

[Pose viewer](../results/diffdock/20260914T145620247011Z/viewer.html) ·
[Raw timing and confidence summary](../results/diffdock/20260914T145620247011Z/summary.json)
