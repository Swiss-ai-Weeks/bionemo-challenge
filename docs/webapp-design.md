# Boltz-2 ranking app — design draft

Status: proposed for discussion. This document does not implement the app.

## 1. Purpose and first version

A small web app for choosing a protein target, selecting multiple small-molecule
ligands, and ranking their predicted affinity with Boltz-2. Users can create
protein variants and compare each ligand against the original sequence and its
variants. Each completed prediction includes a viewable protein–ligand complex.

The primary workflow is:

```text
Choose target -> Add optional mutations -> Select ligands -> Run predictions
                                                               |
                            Compare variants <- View ranking <--+
                                  |
                            Inspect a complex
```

The initial scope is one protein chain and small-molecule ligands expressed as
SMILES. Protein editing supports amino-acid substitutions, including multiple
substitutions in one variant. Insertions, deletions, covalent ligands, multi-chain
targets, model training, and OpenFold3 comparison are deferred.

## 2. Screen structure

Use one workspace with three tabs: **New run**, **Results**, and **History**.
Keep controls visible and avoid a multi-step wizard. On desktop, target and ligand
panels sit side by side. On a narrow screen they stack. Results occupy the full
width, with a complex inspector opening below the selected row.

Visual direction: light background, restrained blue/green accents, readable tables,
small molecule thumbnails, and generous spacing. Color supplements explicit
status labels. Use ordinary buttons and form controls, keyboard-accessible table
selection, and clear focus indicators. The 3D viewer loads only when opened.

All numerical examples in the layouts below are illustrative, not new results.

### New run

```text
+------------------------------------------------------------------------------+
| Boltz affinity lab                                  Model: Ready [green dot] |
| [New run]   Results   History                                                 |
+------------------------------------+-----------------------------------------+
| TARGET                             | LIGANDS                                 |
| [Search name / UniProt ID.......]   | [Search name / ID...................]   |
| [Human Src / P12931          v]     | [Add SMILES] [Import CSV]               |
|                                    |                                         |
| Name: Human Src                    | [ ] Select all filtered (5)             |
| Sequence: 536 amino acids          |     Structure  Name / ID     Reference  |
| [View sequence] [Paste FASTA]       | [x] [diagram]  4521          8.7 nM     |
|                                    | [x] [diagram]  6121          5.1 nM     |
| PROTEIN VARIANTS                   | [x] [diagram]  6122          2.8 nM     |
| [x] Original sequence              | [ ] [diagram]  6123           12 nM     |
| [x] Variant 1: T341I         [Edit] | [ ] [diagram]  6124           25 nM     |
| [+ Add variant]                    |                                         |
|                                    | 3 selected                              |
+------------------------------------+-----------------------------------------+
| Run name: [Src screen.......................................................] |
| 3 ligands x 2 protein sequences = 6 predictions                              |
| Alignment: query sequence only  [Advanced settings v]                        |
|                                                       [Run 6 predictions]   |
+------------------------------------------------------------------------------+
```

- Search the local target/ligand library; do not require external search accounts.
- Seed the library from the existing five-ligand Src example. A background bootstrap
  can fetch the pinned BindingDB archive using the existing downloader. The app
  remains usable for manual input if downloading fails; show a retry action.
- The default ligand list shows compounds associated with the selected target.
  An **All ligands** filter enables cross-target exploration and imported compounds.
- A pasted FASTA creates a new named target. Preserve the exact sequence used.
- Custom ligands need a name and a valid SMILES. CSV import accepts `name,smiles`
  and optional experimental endpoint/value/unit fields, with row-level errors.
- Reference measurements appear only when associated with this target and ligand;
  preserve endpoint, qualifiers, and source. Never silently assign a wild-type
  measurement to a mutant or mix Ki/Kd/IC50 into an experimental ranking.
- Selection persists when filtering. The selection count always includes hidden
  selections; “select all filtered” affects only the current filtered list.
- Show input errors beside the affected control. Disable Run for invalid inputs,
  an empty selection, or an unavailable model, with a specific explanation.
- When no variants exist, Original is selected. When adding a variant, include
  Original by default to support comparison; the user can deselect it explicitly.

### Mutation editor

```text
+------------------------------------------------------------------------------+
| Edit protein variant                                                   [X]   |
| Parent: Human Src / P12931 / 536 aa                                           |
| Variant name: [T341I.......................................................] |
|                                                                              |
| Position*          Original                 Replace with                     |
| [341      ]        T (threonine)             [I - isoleucine              v]  |
| [+ Add substitution]                                                        |
|                                                                              |
| Sequence around position 341                                                |
| Original:  ... [flanking residues] [T] [flanking residues] ...                |
| Variant:   ... [flanking residues] [I] [flanking residues] ...                |
|                                                                              |
| * 1-based position in this exact input sequence, not PDB residue numbering.   |
| Changes: T341I                  Length: 536 -> 536 amino acids                |
|                                                    [Cancel] [Save variant]  |
+------------------------------------------------------------------------------+
```

Offer `T341I` text entry as a shortcut to the same validated substitution table.
Each variant is derived from the immutable original target, not from another
variant. Multiple rows form one combined variant; separate variants are added
explicitly. Reject out-of-range positions, incorrect original residues, duplicate
positions, invalid amino-acid codes, and substitutions that change nothing.

The preview is generated from the actual sequence. Show a numbered full sequence
on demand. Persist the parent sequence, substitutions, and resulting sequence so
saved results cannot change if the target library is later edited.

## 3. Running and ranking

Submit one Boltz-2 request for each selected ligand–protein-sequence pair.
A ligand batch is not a multi-ligand complex. A run snapshots its inputs and
settings before queuing jobs. For one GPU, run jobs sequentially; users can
continue browsing while work runs.

### Live results

```text
+------------------------------------------------------------------------------+
| Src screen                        Running: 4 / 6 finished [Cancel remaining] |
| [====================----------]  Now: Variant T341I / ligand 6121            |
| Protein: [Original sequence v]   [Compare variants]             [Export CSV]  |
+------+-----------+------------------+---------------+-------------+-----------+
| Rank | Ligand    | Pred. IC50 equiv.| Binder prob.  | Confidence  | Status    |
|      |           | nM (lower first) |               | ligand ipTM |           |
+------+-----------+------------------+---------------+-------------+-----------+
|  1   | 6122      |       12.0       |     0.88      |    0.91     | Complete  |
|  2   | 4521      |       29.4       |     0.79      |    0.98     | Complete  |
|  3   | 6121      |       41.0       |     0.75      |    0.89     | Complete  |
+------+-----------+------------------+---------------+-------------+-----------+
| Ranking is provisional while this protein's predictions are incomplete.     |
| Select a row to inspect the predicted complex and measurement reference.     |
+------------------------------------------------------------------------------+
```

Display queued/running/failed jobs in the table with blank scores, never a score
of zero. Results update as individual jobs finish. Keep completed results if any
job fails; provide **Retry failed**. Exact score ties share a rank. Changing the
sort order does not change the recorded affinity rank. Rank within one protein
sequence at a time, not across unrelated sequences.

A compact explanation below the table distinguishes predicted affinity,
probability of binding, and structural confidence. Low confidence is visible in
row details; do not invent a composite score that mixes these quantities.

### Compare original and mutated proteins

```text
+------------------------------------------------------------------------------+
| Variant comparison                         Metric: [pIC50 change vs WT v]   |
| Positive change = predicted stronger binding than the original sequence      |
+-----------+------------------+-----------------+------------------------------+
| Ligand    | Original pIC50   | T341I pIC50     | Change vs original           |
+-----------+------------------+-----------------+------------------------------+
| 4521      |       7.53       |      6.90       | -0.63                        |
| 6121      |       7.39       |      7.10       | -0.29                        |
| 6122      |       7.92       |     Running     | --                           |
+-----------+------------------+-----------------+------------------------------+
| [View affinity values]  [View ranks]  [Export comparison]                    |
+------------------------------------------------------------------------------+
```

For several variants, add one column per variant and allow horizontal scrolling.
Click a cell to inspect that exact ligand–variant complex. Compute a change only
when both predictions succeeded with comparable settings. If Original was not
selected, show absolute values and omit changes. Do not imply a mutation-effect
estimate is experimentally established or provide unsupported error bars.

### Complex inspector

```text
+------------------------------------------------------------------------------+
| Ligand 4521 / Human Src / Original                              [Close]      |
+------------------------------------------------+-----------------------------+
|                                                | Predicted IC50 eq.: 29.4 nM |
|                                                | Predicted pIC50: 7.53       |
|             INTERACTIVE 3D COMPLEX             | Binder probability: 0.79    |
|                                                | Ligand ipTM: 0.98           |
|        Protein cartoon + ligand sticks         |                             |
|                                                | Reference: IC50 8.7 nM      |
|                                                | [Source article]            |
|                                                |                             |
| [Whole complex] [Focus ligand] [Reset view]     | [Download CIF] [Raw JSON]   |
+------------------------------------------------+-----------------------------+
| [Sequence and mutations v] [Prediction settings v] [Interpretation v]        |
+------------------------------------------------------------------------------+
```

Reuse the bundled 3Dmol.js library and existing viewer styling. Highlight mutated
residues on the protein, with a toggle and sequence-position labels. Keep structure
fetching separate from table responses so browsing ranks does not load every CIF.
If the affinity succeeded but a structure is missing or unparsable, retain its
score, show a viewer error, and allow downloading the raw response.

### History

```text
+------------------------------------------------------------------------------+
| History                                             [Search runs...........] |
+-------------------+-------------+----------+-----------+----------------------+
| Run               | Target      | Pairs    | Created   | Status               |
+-------------------+-------------+----------+-----------+----------------------+
| Src screen        | Human Src   |    6     | Today     | Running (4/6)        |
| Src baseline      | Human Src   |    5     | Yesterday | Complete             |
+-------------------+-------------+----------+-----------+----------------------+
| Select a run to reopen results. [Duplicate as new run] copies its selections. |
+------------------------------------------------------------------------------+
```

## 4. Affinity interpretation and defaults

Reuse the current Boltz-2 NIM 1.9.0 adapter and preserve its raw response.
The current container computes its field named `affinity_pic50` with an additional
factor of 1.364. The app must not display that field as dimensionless pIC50.
Use the documented upstream raw-score convention already used by our scripts:

```text
raw score = log10(IC50 in micromolar)
pIC50 = 6 - raw score                      higher = stronger
predicted IC50 equivalent (nM) = 1000 * 10^raw score
rank by raw score ascending                lower = stronger
change in pIC50 = variant pIC50 - original pIC50
```

Version the adapter and conversion rule. Check for finite output scores, preserve
missing values as unavailable, and flag an unsupported model version rather than
silently applying the conversion. Use full precision for ranking and comparisons;
round only displayed values. Experimental measurements are context, never model
inputs. Do not label predictions as Kd or experimentally measured affinity.

Start with the existing demo settings: query-only MSA, 3 recycling steps, 200
sampling steps, 1 structure sample, and 5 affinity diffusion samples with 200
steps. Show these in a collapsed settings panel, with the query-only limitation
visible before running. Apply identical settings to every pair in a run. A later
version can add supplied MSAs and repeat runs to estimate ranking stability.

## 5. Lightweight architecture

Proposed stack: Python + FastAPI, server-rendered Jinja templates, plain CSS, and
small vanilla JavaScript modules. Use SQLite for metadata and a filesystem volume
for requests, responses, and CIFs. Use RDKit only on the backend for SMILES
validation/canonicalization and SVG ligand thumbnails. No frontend build step,
Node runtime, Redis, Celery, or separate database container is needed.

```text
Browser
  | HTML/CSS + small JS       | lazy-load bundled 3Dmol.js
  | poll JSON every 2s while a run is active
  v
+---------------------------- app container -------------------------------+
| FastAPI + Jinja                                                          |
|   | validate input / create run / retrieve ranking / serve artifacts       |
|   +-------------------------> SQLite                                     |
|                                  ^                                       |
| One background worker thread ----+ claim queued job                      |
|   |                                                                      |
|   | one request per ligand + protein variant                             |
|   +---------------------> http://boltz2:8000/.../predict                   |
|   |                                                                      |
|   +----> save request / response / CIF on persistent volume               |
+--------------------------------------------------------------------------+
                                      |
                      +---------------v---------------+
                      | Boltz-2 NIM container          |
                      | One GPU + persistent cache     |
                      +-------------------------------+
```

Run exactly one application process with one worker thread in v1. The worker uses
its own SQLite connection and short transactions; no database transaction stays
open during GPU inference. Polling avoids WebSocket infrastructure. This is a
local research workspace, with no user accounts or collaboration features in v1.
Bind the published app port to localhost by default; remote users can use SSH
forwarding. Supporting public/shared deployment would require a separate auth
and access-control design.

### Persistence and job lifecycle

```text
Target --< Variant
Ligand
Run --< Job >-- immutable target/variant + ligand snapshot
          |
          +--> request.json / response.json / complex.cif / error details

Job: queued -> running -> succeeded
                       -> failed
                       -> interrupted (application restarted)
     queued -> cancelled

Run: queued / running / complete / completed with errors / cancelled
```

Store the parent sequence, mutated sequence, ligand SMILES, source identifiers,
settings, model version/digest, timestamps, raw score, conversion version, and
artifact paths. Editing a library item never changes an existing run.

Mark in-flight jobs interrupted on app restart; do not automatically resubmit an
uncertain GPU request. Pending jobs remain durable, and resume once the worker
has confirmed the backend is idle. If backend completion is uncertain, require
an explicit retry rather than overlap a second request. Surface timeouts and
retry attempts in job history. Retry failed/interrupted jobs creates new attempts
without erasing previous diagnostics or successful results.

**Cancel remaining** cancels queued jobs and lets the active GPU request finish;
label this behavior directly. Stop requesting new jobs when cancellation is set.
Refreshing or closing the browser does not cancel the run. Do not cache predictions
across runs in v1: repeating a run means fresh predictions.

### API outline

```text
GET    /api/status                       model readiness and worker state
GET    /api/targets                      local search
POST   /api/targets                      add sequence / FASTA
POST   /api/targets/{id}/variants        validate and save substitutions
GET    /api/ligands                      filter/search local library
POST   /api/ligands                      add SMILES
POST   /api/ligands/import               import CSV, report invalid rows
POST   /api/runs                         snapshot selections; enqueue jobs
GET    /api/runs                         history
GET    /api/runs/{id}                    status, ranked scores, comparisons
POST   /api/runs/{id}/cancel              cancel queued work
POST   /api/runs/{id}/retry               retry selected failed attempts
GET    /api/jobs/{id}/structure           CIF for viewer
GET    /api/jobs/{id}/artifacts/{kind}    allowlisted download types
GET    /api/runs/{id}/export              ranking/comparison CSV
```

Serve artifacts by database identifiers and allowlisted kinds, never arbitrary
filesystem paths. Validate requests on the server, bound sequence/batch sizes,
and render imported names as text. Keep the NGC key inside the backend deployment;
it never goes to the browser, exported runs, or application logs.

## 6. Docker Compose experience

With the existing `.env` key, registry login, and GPU Docker setup in place:

```text
docker compose up --build
                 |
                 +--> cache-init: prepare volume ownership, then exit
                 +--> boltz2: download/cache model, warm up, become ready
                 +--> app: start immediately, show model initialization status

Open http://localhost:8080
```

The steady-state command is simply `docker compose up`; Compose builds the app
image on first use. Use `--build` after changing its image contents.

Add an `app` service with a local Dockerfile and publish `127.0.0.1:8080:8080`.
Use the Compose service name `boltz2` for inference, not host-loopback addresses.
Make its GPU index configurable with a default of GPU 0, so the app needs only
one suitable GPU. OpenFold3 stays behind the existing optional profile and is
not required by this app; on a two-GPU host it must use a different GPU when enabled.

The app starts without waiting for model health, shows **Model starting**, and
disables submission until ready. Do not make a lengthy model warmup look like a
broken website. Persist application state and artifacts in an `app-data` volume;
retain the existing Boltz model/output volumes. No manual Python commands are
required for normal operation. Dataset bootstrap runs once in the background,
records its status, and can be retried from the UI.

Keep `.env`, downloaded datasets, generated inputs/results, and SQLite files out
of Git. Include code, templates, static assets, dependency pins, and documentation.

## 7. Delivery slices and review points

1. **Workspace and library:** Compose app service, target and ligand selectors,
   mutation validation, example bootstrap, and model readiness display.
2. **Prediction runs:** persistent queue, one pair per request, progress, failure
   recovery, ranking, and CSV export.
3. **Inspection and comparison:** interactive complexes, mutation highlighting,
   original/variant comparison, and history.

Acceptance checks should cover mutation indexing and sequence integrity; correct
ligand-by-variant job counts; affinity conversion and ties; persistence after a
restart; one failing job without losing other results; and a real small batch
against the local NIM. Browser checks should cover selection, mutation errors,
progress, sorting, and viewing the selected ligand–variant complex.

Proposed defaults for our next iteration:

- One-page setup instead of a wizard.
- Original plus named substitution variants; no indels initially.
- Ranked table first, variant-comparison table second, complex inspector on demand.
- A single app process with SQLite and one GPU worker.
- Local deployment with one required GPU and `docker compose up`.
