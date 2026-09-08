# Assay conversion repair

The compatible `normalize(value, unit, relation='=', endpoint='IC50')` API now
retains endpoint identity, reverses concentration inequalities under the decreasing
negative logarithm, handles all six required units and refuses invalid inputs.
The defining quantity is `p_activity = -log10(concentration / 1 M)`. A censored
threshold remains a censored threshold, not a point observation. This helper does
not establish assay comparability or whether a row is eligible for model fitting.

The interface contract is `inputs/legacy/README.md`; the adjacent caller inspected
was `inputs/legacy/vendor_model.py`. The original source is unchanged. Successful
return keys and defaults are preserved. Unit spelling remains case-sensitive;
Unicode ≤/≥ comparisons canonicalize to ASCII <=/>=. `IC50`, `Ki`, and `Kd` are
the supported, distinct endpoints. Numeric strings are supported. Booleans are
refused explicitly rather than silently interpreted as 0/1 concentrations. A
concentration and its nM conversion must be finite and positive in binary64.
Extreme values that overflow the nM result are rejected. No concentration floor
is invented. Arbitrary nonnumeric objects and unsupported metadata are invalid;
this is not a general dimensional-analysis parser.

## Frozen references and observed results

The author checks use analytic powers of ten and explicit reversed relations,
not a round trip through the candidate. Before running the candidate, expected
semantics and tolerances were fixed in `analysis/test_assay_math.py`: absolute
p-activity tolerance 1e-12 and nM relative tolerance 1e-12. These tolerances cover
binary64 evaluation in the selected simple examples; categorical equality is
exact. Equivalent unit cases all represent 1 µM (1000 nM, p=6); 10 nM has p=8;
100 nM has p=7. Tiny and large positive finite cases test removal of the original
floor without arbitrary clipping. These are author tests and are not blinded.

- Immutable baseline SHA256: `5b5d74c6577f169bd6508222e3dbee69b7202791e87dc5c7d93d2d03b1b16f92`.
- Candidate SHA256, fixed before independent scoring: `7b60a73af904b265beaad53ec7ef03450847deff9a95ccfa7bfb6d72ccaaf941`.
- Test-source SHA256: `0be8a08f307867da2ab8d71482157762acb8d7d3e59351b9e8d1290473bc31aa`.
- Baseline: 39 collected, 10 passed, 29 failed, zero skips/xfails; command exit 1.
- Candidate revision 1: 39 collected, 39 passed, zero failures/skips/xfails; exit 0.
- Exact case outputs and actual imported module/interpreter paths are retained in
  `analysis/checks/assay_math_baseline.json` and `assay_math_candidate.json`.

Baseline failures include missing micro-symbol/mM conversions, unchanged
inequality directions, Ki/Kd relabeling, a floor applied to valid tiny
concentrations, acceptance of zero/negative/nonfinite values and unknown metadata,
and uncaught unsupported input exceptions. The ten passing baseline cases are
retained as positive controls. No failed candidate iteration was discarded;
candidate revision 1 passed its author tests on the first execution. Root's
separately authored oracle passed 12/12 checks against the frozen candidate:
`analysis/checks/assay_independent_oracle.json`, with replay source
`analysis/independent_assay_oracle.py`. It supplements author cases using exact
molarity-power references and a separately expressed ln(2)/ln(10) reference,
reversed bounds, endpoint retention, and seven rejection cases. The oracle
shares Python's math primitives and the supplied scientific definition, but its
authorship is independent. Both source and contract were visible; this is
non-blind independent evidence, not hidden evaluation.

## Replay and reconstruction

Use the pinned common runtime in `admin/ENVIRONMENT.json` and an authorized CPU
lease. Original execution used CPU3, one thread, no GPU, and only two short Python
test processes; no background process was started. Each command completed in
under one second. Memory peak was not measured; the standard-library-only test
is bounded to 39 scalar cases. Installed environments were not modified.

```sh
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 taskset -c 3 /home/zbos/.local/share/cedar-runtime/20260908/cedar-python analysis/test_assay_math.py --module inputs/legacy/assay_math.py --output NEW_OUTPUT/assay_math_baseline.json
env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 taskset -c 3 /home/zbos/.local/share/cedar-runtime/20260908/cedar-python analysis/test_assay_math.py --module analysis/assay_math.py --output NEW_OUTPUT/assay_math_candidate.json
patch --output=NEW_OUTPUT/assay_math_reconstructed.py inputs/legacy/assay_math.py analysis/assay_math.patch
cmp analysis/assay_math.py NEW_OUTPUT/assay_math_reconstructed.py
```

Create `NEW_OUTPUT` before replay and select a currently authorized CPU. Baseline
exit 1 is expected and must not stop the candidate check. The delivered unified
patch was generated with `diff -u --label inputs/legacy/assay_math.py --label
analysis/assay_math.py inputs/legacy/assay_math.py analysis/assay_math.py` (diff exit
1 denotes differences), reconstructed into `/tmp/cedar-assay-math-reconstructed.py`
with `patch`, then compared byte-for-byte using `cmp` (exit 0). Its SHA256 matched
the frozen candidate. No unrestricted or hidden-test validity is claimed.

## Actual guidance use

Read `.agents/skills/scientific-software/SKILL.md`, its
`references/scientific-invariants.md`, and
`.aidd/expertise/knowledge/methods/scientific-software.json`. These informed the
physical-quantity contract, immutable baseline evidence, analytic references,
boundary/negative controls, candidate identity and reconstructable patch, and
request for independent oracle results. Also read the expertise
`docs/GETTING_STARTED.md` entry path. No linked external sources, prior trials or
other repositories were opened. AIDD expertise use is recorded; AIDD captured
execution was not configured or claimed for this repair.
