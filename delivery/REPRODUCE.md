# Reproduce the checked delivery

The immutable input manifest SHA-256 is `77934c089aa941e191e08c2459d12d6c522596ba4614f81c7fb7850aba6a7094`. All 30 member file hashes are verified before curation and final validation. Inputs must remain byte-identical. The code is original analysis plus a compatible repair of the supplied simulated helper; its reconstruction patch and attribution are preserved. Source data terms are in `inputs/project/PROVENANCE.md` and supplied OpenFF license data. Scientific internet is closed; no network, model download or package installation is needed.

The exercised runtime is Python 3.12.14, RDKit 2026.03.6, NumPy 2.5.2, SciPy 1.18.0, scikit-learn 1.9.0, pandas 3.0.5, matplotlib 3.11.1, PyMOL 3.2.0a. The installed common launcher is `/home/zbos/.local/share/cedar-runtime/20260908/cedar-python`; its base interpreter is `/home/zbos/.cache/molmonth-tools/gromacs-env/bin/python`. RDKit/PyMOL use the common overlay and other listed packages use that base environment. Preflight import checks and structural evidence retain observed version/path details. Pillow and the system DejaVu Sans font support static image composition; the molecular render receipt records the font path/hash. No installed environment or inventory source file is included in Git. On another host, supply a Python executable/launcher with these dependencies via `CEDAR_PYTHON`; changes require revalidation rather than presumed numerical equivalence.

From the repository root, with an approved one-core / 1 GiB resource lease and its allowed CPU ID, run (replace `2` if a later grant assigns another CPU):

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export CEDAR_PYTHON=/home/zbos/.local/share/cedar-runtime/20260908/cedar-python
export MPLCONFIGDIR=/tmp/cedar-matplotlib-replay
taskset -c 2 "$CEDAR_PYTHON" analysis/reproduce.py --out /tmp/cedar-replay-new
```

The output directory **must not already exist**. The driver copies approved immutable inputs, analysis source/notes and report documents, then executes against those copied inputs. Original outputs are preserved. It runs curation, original vendor reproduction, common-fold model audit, actual legacy split-overlap check, experiment accounting, structure audit, inventory/source tautomer diagnostic, helper tests/independent oracle, static figures and the final validator. It does not create a Git repository, publish, contact another task or access scientific sources outside the bundle. Each child inherits the CPU affinity and one-thread limits. The expected duration is a few minutes with the frozen runtime, dominated by two static PyMOL views; no GPU is used. `--skip-molecular-render` supports a numerical-only replay but is explicitly incomplete for figure regeneration (original images are copied and labeled as such).

Expected numerical outputs:

- Curation: 266 activity rows; 122 eligible, 138 context, two duplicate and four quarantined. 258 source comparisons all match. 28 request predictions: 18 observed, six abstained, two duplicate and two quarantined.
- Original vendor: 262 rows, R²0.9601857410, MAE0.0668989289. Actual random split shares one source activity ID and nine chemical graphs. This is intentionally invalid prospective evidence; baseline defects remain reproduced.
- Common scaffold folds: 122 rows / 117 chemical graphs / 26 scaffold groups. Descriptor RF MAE0.5390382641 and R²−0.6025680057; retrospective-index RF MAE0.0422035043 and R²0.9842803150. No fit is deployed for cellular inference.
- Cycle: $17,800 total; eight unique compounds, six advances/two diagnostics; 25 combinations/26 charged experiments; latest completion day7; every material check passes.
- Original export fixed-frame RMSDs: reference0.000000Å, priority8.774964387392Å. Priority contacts unusable. Strict graph/state maps and independent scalar/frame controls pass. CED23 exact state differs, canonical-tautomer parent matches.
- Candidate helper: 39/39 author checks and12/12 independently authored analytic checks. Immutable baseline:10/39 passes and29 preserved failures. The patch reconstructs the candidate byte-for-byte.
- `analysis/validation.json` and `analysis/replay_receipt.json` report actual validation and child exits. Required numeric comparisons use independent Decimal molarity/log calculations and per-vial cost/material/time accounting. A schema/number pass does not validate future biology.

The driver preserves baseline failure logs as expected failures; any unexpected child failure stops the replay and leaves partial output with a failure receipt. Model seeds and folds are fixed, so numerical tables should reproduce under these versions. Absolute path metadata, timings, rendered-image binary encodings and saved PyMOL session internals may differ without a scientific quantity change; do not claim byte-identical environments or figures. Input and tested helper hashes must match. Native captured execution is absent; the files are ordinary reproducible scientific work.

For a narrow helper reconstruction, apply the patch in a fresh directory containing the immutable base (see `analysis/assay_repair_notes.md` for the exercised exact command), then compare SHA-256 and run the tests. The original source remains under `inputs/legacy/`. The legacy launcher does not support Python stdin `-`; use the retained script files or `-c`, not the failed stdin pattern documented in progress notes.
