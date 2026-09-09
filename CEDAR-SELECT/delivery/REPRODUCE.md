# Reproduce this delivery

Use the final Git commit in the completion message. Source inputs and final scientific outputs are bound by `MANIFEST.sha256`, whose paths are relative to repository root. Its own digest is reported separately to avoid self-reference. No installed environment, private guidance, native trace or cache is required for replay.

Prerequisites: pinned launcher `/home/zbos/.local/share/cedar-runtime/20260908/cedar-python`, Python3.12.14, RDKit2026.03.6, NumPy2.5.2, SciPy1.18.0, scikit-learn1.9.0, pandas3.0.5, matplotlib3.11.1, Biopython1.88, PyMOL3.2.0a. The current runner uses only a subset of these installed packages. No packages were installed or upgraded. Exact import/runtime records are in operations/preflight.md and per-job JSON. On another machine provide compatible ordinary packages and set CEDAR_PYTHON to the interpreter/launcher path; cross-platform graphical byte identity is not claimed.

From a clean checkout's root, verify the frozen package:

```sh
sha256sum -c CEDAR-SELECT/delivery/MANIFEST.sha256
```

To rerun without overwriting reference artifacts:

```sh
python3 CEDAR-SELECT/delivery/replay.py /tmp/cedar-replay-report
```

This copies only CEDAR-SELECT and the predecessor inventory needed for identity matching into a fresh temporary directory, deletes generated accepted outputs there, runs nine sequential jobs with one OMP/BLAS/MKL thread and600-second per-job timeout, compares scientific JSON excluding runtime fields and all seven model CSV files byte-for-byte, verifies source bytes, and removes the temporary computation directory. Logs/report remain at the requested output path. A source derivative of the supplied PDF/validation text is already present; extraction used pdftotext, and Table2 was visually checked. Rendering is regenerated but scientific conclusions do not depend on identical PNG compression.

Executed final replay: all nine commands exited0; nine scientific JSON comparisons and seven CSV comparisons passed in8.28s total wall. No scientific outputs in the reference project were overwritten. Evidence: [replay_report.json](../operations/replay/replay_report.json). Model numerical checks and independent source-level audit must also pass. This is same-runtime reproducibility and nonblind independent numerical validation, not unseen-data model validation.

Entry points and outputs:

- model-work/run.py: curation, source/group partition audit, fixed models/CV/development, null and cluster stress; outputs model-work/results.
- coverage-work/analyze.py: analogue transformations/crosswalk, typed retrieval audit, obligation/material/timing accounting.
- design-work/designs.py: exact isomeric structure and source matching, static descriptors; reads inputs/data/inventory.csv solely as predecessor chemistry.
- model-work/score_designs.py: training-only diagnostic outputs and explicit same-core/similarity checks; failed model gate forbids numerical advancement.
- coverage-work/model_check.py: independent audit without importing root model pipeline; retained limitations include tiny cluster support and descriptive uncertainty.
- structure-work/analyze.py, sensitivity.py and render.py: accepted sequence-based contact/state analyses and molecular view. Do not substitute analyze_v1_same_number.py or its rejected output.
- model-work/figures.py: analogue tradeoff and two-dimensional chemistry views.

The default scripts write to their own subdirectories, so execute them directly only in a disposable replay checkout, not in a frozen delivery. `replay.py` is the safe default. It does not reproduce the rejected exploratory attempts or live literature searches, whose failures/source dates remain recorded. Scientific model predictions here are diagnostic and fail the advancement gate; replay PASS never changes that conclusion.

For a portable source bundle, archive the final Git commit's CEDAR-SELECT tree plus inputs/data/inventory.csv, inputs/MANIFEST.json and admin/ENVIRONMENT.json. Do not include `.aidd`, `.agents`, `.codex`, `.git`, environments, caches, credentials or native traces. Existing public repository access is authorized; no other publication was performed.
