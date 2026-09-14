# Example and interpretation

The first BindingDB Articles row supplies human Src (P12931, 536 residues) and
ligand 4521 (SMILES). The recorded experimental IC50 is 8.7 nM. This demo uses a
query-only alignment and the full database sequence; the assay construct has
not been checked. Use homologous alignments and verify the construct before
interpreting results quantitatively. No affinity labels are supplied to either model.

OpenFold3 predicts structure and confidence. Boltz-2 independently predicts a
complex and affinity; it does not score the OpenFold3 output in this workflow.
In NIM 1.9.0, the field called `affinity_pic50` is calculated as
`(6 - affinity_pred_value) * 1.364` in the container code, so it is not
dimensionless pIC50. The summary instead follows upstream Boltz's convention:
raw score = log10(IC50 in micromolar), pIC50 = `6 - raw`, and IC50 in nM =
`1000 * 10 ** raw`. Preserve raw scores and do not mistake structural confidence
or binding probability for affinity. See `results/boltz2/affinity_summary.json`.

The version-specific `affinity_pic50` observation above is recorded from the original local run; recheck it when upgrading the NIM image.

Request documentation: [OpenFold3](https://docs.nvidia.com/nim/bionemo/openfold3/latest/example-requests.html) and [Boltz-2](https://docs.nvidia.com/nim/bionemo/boltz2/latest/inference.html).
