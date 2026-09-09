# Independent model artifact validation

Separate script `model_check.py` reads the supplied source JSON and root artifacts without importing or executing `model-work/run.py`. No models were fitted or tuned. The checked hashes and actual resource use are in `model_validation.json`. Supplied-data checks pass; the model-led advancement rejection is supported by the recalculated primary metrics.

All 120 source rows have unique IDs and valid target/label state: 105 exact labels, 15 CDK1 pIC50≤5 bounds across 60 compounds. All paired rows agree on isomeric graph, partition, group, study and series. All source-group, nonisomeric-graph and fingerprint collisions are kept together within 56 merged groups; no group crosses training/development. All training fold assignments cover eligible training compounds once and keep related compounds in one fold. Primary populations are CDK2 49 train/11 development and CDK1/ratio 36 train/9 development. Bounds are absent from exact regression populations and retained with correct inequality reversal for concentration and ratio. Hardcoded bound=5 in root is correct for these supplied bytes but should not be generalized to a changed input schema.

Every original-unit analogue transform and context crosswalk was independently checked: the same 16 labels recur under different source identifiers; they were not appended. For the other context labels the source JSON supplies pIC50 rather than original concentration, so this audit verifies consistent reuse and reverse transformation, not original publication transcription. All 60 concentration/ratio transformations, 15 bound residuals and 16 crosswalk records agree.

All nine endpoint/model development MAE, RMSE, bias, rho and counts agree to 1e-10. Every mean and 3NN prediction and every nearest-training Tanimoto value were independently reconstructed using source-order tie breaking and RDKit similarity. Ridge development metrics were independently verified from stored predictions, not refitted. Model selection agrees with training CV metrics and declared tie order; training OOF predictions were not exported, so this audit does not independently rederive CV fits or descriptive error quantiles.

| Endpoint | Selected | Development MAE | Mean baseline MAE | rho | Gate |
|---|---|---:|---:|---:|---|
|CDK2 pIC50|3NN|0.299395|0.241551|0.590909|Fail: error worse than mean|
|CDK1 pIC50|mean|0.252649|0.252649|undefined|Fail|
|log10 selectivity|mean|0.244966|0.244966|undefined|Fail|

The consequential stress-test limitation is endpoint-specific. Structural clusters have sizes **57,1,1,1**. CDK2 cluster predictions cover all 60 compounds, but 57 predictions train on only three compounds. CDK1 and log-selectivity cluster predictions cover **only one of 45 exact-label compounds each**: the dominant cluster is skipped because too few eligible training compounds remain. Their small MAEs cannot be presented as successful scaffold extrapolation. Retain the one-point values as descriptive outputs with 44 unscored exact compounds and an explicit stress-test abstention. No unseen study/scaffold qualification is earned.

The root protocol's domain rule (nearest similarity≥0.65 and same spiro core) is not implemented as candidate-domain acceptance in this evaluation script. That is not a defect in these fixed development error calculations, but later candidate outputs need a separate actual core/domain check. Model disagreement and 90th-percentile grouped OOF absolute errors are descriptive, not calibrated uncertainty; using selected-model OOF errors adds selection optimism. Protocol claims correctly remain nonblind, same-public-study development. Model failure does not itself prove a chemistry direction should stop.

Three independent audit negative controls were rejected: a reversed bound inequality, a bound decorated with an exact pIC50 value, and a µM concentration treated as nM. These exercise this audit's invariants, not mutation runs through the root pipeline. No live retrieval, hidden labels, external computation, installs or environmental changes. One numerical thread, peak RSS approximately 130 MiB, under one second CPU. All jobs stopped; no subdelegates.
