# Coverage evidence and actionable allocation

The supplied public analogues favor a measured fluoroethoxy comparator (S3, patent example 39) over paying now for unverified S4 procurement. They do **not** establish that successive fluorination improves CDK2/CDK1 discrimination. Recommend proposing O2 and O1, maintaining W2, and retaining $2,600 of the $10,000 component. All logistics are simulated; no expenditure, contact, assay or procurement was executed.

## Assay scope, overlap and local evidence

`analogue_reconciled.json` preserves all eight full CX identifiers, isomeric structures, original observation IDs and original values. It verifies all 16 source pIC50 transforms and all four catalogue bindings by stereochemistry-preserving canonical SMILES. All eight compounds and all 16 target labels match the supplied public context exactly (within 1e-10 pIC50), under different HC/HR identifiers. These are overlapping source observations, **zero additional independent measurements**. There are eight compounds, seven supplied leakage groups, one study, and no reported independent repeat count. Ex69/70 remain separate stereoisomers despite sharing a leakage group. No labels are inferred from simulated work completions or retrieval entries.

Public evidence source: supplied attributed extract of US20240360137A1, with patent example/page/source-panel pointers retained. No external retrieval was required for this bounded supplied-record analysis. Paired IC50 context is human-attributed CDK2/E1 versus CDK1/B1, 1 mM ATP, TR-FRET RB1 Ser780, 25°C/60 min, 100 versus 10 pM enzyme; constructs and source processing must be preserved. Ratio means CDK1 IC50/CDK2 IC50 in these conditions, not physiological selectivity. Configuration and sequence attribution are source-assigned, not newly established.

| Example / slot | CDK2 nM | CDK1 nM | Ratio |
|---|---:|---:|---:|
|12 / S1, methoxy|20|2240|112|
|23 / S2, ethoxy|15|2430|162|
|39 / S3, fluoroethoxy|7|1200|171.43|
|41 / S4, difluoroethoxy|5|730|146|
|56, trifluoroethoxy|10|800|80|
|69, difluorocyclopropoxy stereo A|6|800|133.33|
|70, difluorocyclopropoxy stereo B|6|820|136.67|
|77, hydroxy|121|5820|48.10|

The local attachment environment is pyrazolyl O substitution on the same bicyclic lactam/spirocyclopropyl core with unchanged source-assigned hydroxycyclohexyl stereochemistry. `local_pairs.json` retains six explicit before→after comparisons and all four observation IDs per comparison. Delta pIC50 is after minus before; positive means stronger inhibition. O-methylation 77→12 gives +0.782 CDK2 and +0.415 CDK1 (ratio 48→112), while methoxy→ethoxy 12→23 gives +0.125 and −0.035 (112→162). Each specific transformation has one pair, not replicated support.

For terminal ethoxy fluorination in this narrow environment, 23→39 supports a small ratio increase (162→171); 39→41 and 41→56 oppose a general monotonic-selectivity claim (171→146→80). Thus 1 supporting and 2 opposing pair directions across four unique compounds; shared endpoints make them correlated, and progressive fluorination changes the immediate chemical environment. They must not be pooled as three independent trials of an identical edit. On CDK2 alone, the first two support stronger inhibition (+0.331,+0.146), but the third opposes it (−0.301). Neither direction has a uncertainty-qualified numerical gain. The 69→70 stereo pair has equal rounded CDK2 potency and only −0.0107 CDK1 delta: one unresolved pair for a meaningful stereo advantage, not proof of equivalence. No incompatible observations were included in these numerical comparisons.

Use S2/S3 as paired prospective benchmarks; do not spend a scarce early slot chasing a 7→5 nM difference while CDK1 inhibition also strengthens. Reversal: reproducible, bracketed paired-context data showing S4 has materially better discrimination than S3, or verified exposure/solubility benefit sufficient to improve a subsequent cellular window, would justify reconsidering S4. Absent such evidence, fluorination is a potency/physicochemical hypothesis, not a demonstrated selectivity solution.

## As-of budget, material and schedule

By obligation identity, F1/F2 are one paid OB-17 ($1,800); F3/F4 one reserved OB-23 ($2,200). $6,000 remains before new actions. The entire $10,000 including this $4,000 is within the $120,000 tranche; the rest of the tranche is $110,000, not $120,000 in addition. No handling charge is added to all-in option totals.

By vial identity, V1 has 5 mg free, V2 has 6 mg total including 2 mg already reserved (4 mg free; I2/I3 are duplicate exports), V3 has 4 mg free. Consumed amounts are already excluded and must not be subtracted again. M3 confirms simulated received V3 stock and is already funded in OB-17. M4 attests one historical simulated S1 synthesis, not additional stock, current potency or S4 feasibility.

W1 is valid completed S1 CDK2 but its number is unavailable. W2 has only CDK2 valid completed; its CDK1 endpoint is pending day 3 with 2 mg V2 reserved. W3 is failed CDK1 QC, no eligible potency estimate. W4/O1 is an explicitly intentional independent-preparation repeat, not an accidental duplicate and not already authorized.

| Proposed option | Incremental | Material | Timing | Rationale |
|---|---:|---|---|---|
|O2, S3 CDK1 repeat with QC|$2,000|2 mg V3|lane 2, day 0–3|Repair failed selectivity endpoint for useful available comparator|
|O1, S1 CDK2 independent repeat|$1,400|1 mg V1|lane 2, day 3–5|Challenge preparation/day reproducibility of baseline|
|W2 CDK1 pending|$0 new (reserved $2,200 obligation)|2 mg V2 reserved|conservatively lane 1 through day 3|Preserve existing authorization; obtain valid pending endpoint|

This schedule never exceeds two concurrent jobs and finishes by day 5, within day 8. It assumes the pending work occupies one lane until its recorded due date; no hidden cancellation or optimistic overlap is required. Residual free stock is V1 4 mg, V2 4 mg, V3 2 mg. Total earmark spend/reservations/proposed commitments = $7,400; remaining $2,600 is held, not silently allocated elsewhere. Pending delays still permit O1 through day 5 on lane 2. No option above creates a new valid assay result until QC passes.

O3 costs $1,600 and fits stock but duplicates S2's pending CDK1 endpoint without an independent-repeat purpose. Defer. If W2 fails QC on day 3, O3 may be newly proposed as an explicitly identified recovery assay (2 mg V2; day 3–6 on lane 1, after confirming QC correction and scope); it would leave $1,000. Do not cancel OB-23 or recharge its existing work. This conditional spend is not included as an authorization now. O2 alone is a lean alternative retaining $4,000; O1 earns $1,400 because W1 supplies no accessible numeric result and preparation/day repeatability remains unresolved. Before interpreting O1, obtain W1's numerical result and protocol/QC record; without that comparison it is one standalone assay, not a measured repeat spread.

O4 is not budget- or deadline-qualified: price, lead time and vial are null. M1 is a two-step retrosynthetic proposal without procedure/yield; M2 is an unverified listing without quantity, price or lead time. Neither is demonstrated S4 synthesis nor purchasable stock. Require identity/purity, >=2 mg available quantity, all-in quote within then-uncommitted $2,600, and arrival plus panel duration by day 8 with lane capacity. Confirm those facts before proposing commitment; scientific value must also beat retained reserve. No cost or completion date is invented.

O2 passing QC alone does not establish S3's new paired selectivity because no new S3 CDK2 endpoint is supplied. Retrieve a documented comparable paired CDK2 result or propose that endpoint under a costed later protocol before calculating a new ratio. W2 needs its missing CDK2 numerical/QC record alongside its CDK1 completion. Progress preference is CDK2 ≤100 nM and ratio ≥10 in documented paired conditions; QC failure, unbracketed curves or uncertain context blocks advancement rather than being called inactivity. A reproducible failure of these biochemical criteria reverses comparator advancement; cellular advancement additionally needs target engagement and a CDK1/normal-cell window, outside what these options establish.

## Retrieval and typed-entity audit

`retrieval_audit.json` accounts for all nine envelopes, their intent, final executed query, status, page totals, typed entities and record-level restrictions. R1 is an eligible simulated index hint only (no numerical measurement). R2 executes unsupported `target` instead of supported `target_accession`, and returns P06493 rather than P24941: repair query, exclude this row from CDK2 evidence. R3 changes compound to assay namespace; token 771 equality is not identity and no assay→compound crosswalk exists. R4 is partial because advertised page 2 is absent. R5 is a completed supported empty lookup for compound 999, not biological inactivity. R6 is HTTP 503 failure, not empty. R7 is valid empty form mechanism lookup; documented indexing/relation permits bounded R8 parent lookup. R8's parent annotation is not a child assay label, proof of actual salt identity or potency. R9's percent inhibition at 10 µM fails IC50 eligibility despite correct query intent; no IC50 conversion is justified. All available returned records were checked, not sampled. Fictional service contracts are the authority for this exercise; no live provider behavior or search capability was tested.

## Replay and resources

Run `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 /home/zbos/.local/share/cedar-runtime/20260908/cedar-python CEDAR-SELECT/coverage-work/analyze.py` from project root. The script writes only coverage-work, checks transforms and binding identity, deduplicates financial/inventory records with conflict assertions, and retains exact source IDs. Resource log records actual CPU/wall/RSS and input hashes. No dependency/environment changes, live retrievals or subprocess jobs. All work here has exited; no jobs/delegates remain for this subtask. Guidance used: potency-modeling local analogue, assay/retrieval integrity and shared scientific selection procedures, informing context-specific pair interpretation, typed joins and conservative commitment accounting. Private guidance content is not copied into outputs.
