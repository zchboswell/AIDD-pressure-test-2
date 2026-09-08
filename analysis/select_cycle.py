"""Explicit decision manifest and auditable hard-constraint accounting."""
import argparse,collections,json,pathlib
from curate import readcsv,writecsv
ROOT=pathlib.Path(__file__).resolve().parents[1]
ADVANCE=['CED-02','CED-03','CED-08','CED-13','CED-17','CED-22'];DIAGNOSTIC=['CED-23','CED-21']
REASONS={
'CED-01':'Reserve: TE65 nM coexists with CDK1 90 nM and half-life14min; CED21 represents this unfavorable profile.',
'CED-02':'Advance conditionally: TE130 nM, solubility18uM, half-life22min; determine matched biochemical CDK1 window; stability threshold proximity remains uncertain.',
'CED-03':'Advance conditionally from VIAL-03 only: favorable solubility35uM/half-life38min; TE and CDK1 missing. IDENTITY must precede biology because STOCK-28 conflicts; it cannot reassign historical data.',
'CED-04':'Reserve for a later cycle: arrival10business days exceeds deadline; TE850nM also weak.',
'CED-05':'Reserve: purity72percent requires IDENTITY before any new biology; TE35nM does not offset CDK1 45nM/half-life8min. No assays allocated; purity remains unresolved.',
'CED-06':'Reserve: low solubility3uM, half-life14min and CDK1 90nM; CED21 diagnostic covers unfavorable profile.',
'CED-07':'Reserve: same TE130/CDK1 400/half-life22 as selected02/17/22 but poorer8uM solubility; prioritize more tractable peers.',
'CED-08':'Advance exploratory: CDK1 1800nM, half-life38min, solubility18uM; cellular engagement missing. Two independently prepared ENGAGE runs resolve a consequential unknown.',
'CED-09':'Reserve: TE850nM misses preference despite good55min stability; missingCDK1 alone does not justify slot over selected missing-TE candidates.',
'CED-10':'Reserve: TE35nM with CDK1 45nM and half-life8min; unfavorable comparator already represented.',
'CED-11':'Reserve: TE missing but half-life14min/CDK1 90nM less promising than CED08; do not infer activity from missingness.',
'CED-12':'Reserve for later cycle: arrival10business days exceeds deadline, and solubility3uM is poor.',
'CED-13':'Advance as explicit tradeoff experiment: TE320nM is weaker than preference; CDK1 1800nM/half-life38min merit matched-panel test. No potency pass unless new data meet gate; solubility8uM remains concern.',
'CED-14':'Reserve: TE850nM with missingCDK1; prioritize smaller cellular gap in13 rather than stable-but-weak candidate.',
'CED-15':'Reserve: TE35nM with CDK1 45nM/half-life8min; favorable potency alone is insufficient.',
'CED-16':'Reserve: missingTE but half-life14min/CDK1 90nM; lower information value than08.',
'CED-17':'Advance conditionally: TE130nM, solubility120uM, half-life22min; matched-panel comparator to02/22; repeatability remains unresolved.',
'CED-18':'Reserve: sameTE320/CDK11800/half-life38 as13 but solubility3uM vs8uM; prioritize13.',
'CED-19':'Reserve: CDK16000nM/half-life55min but TE850nM;13 tests a smaller cellular potency compromise.',
'CED-20':'Reserve: TE35nM with CDK145nM/half-life8min; unfavorable profile already represented.',
'CED-21':'Diagnostic unfavorable-profile comparator: TE65nM with CDK190nM/half-life14min. New matched panel can confirm or overturn presumed poor window; not a proven nonselective standard.',
'CED-22':'Advance conditionally: TE130nM, solubility70uM, half-life22min; define matched window; borderline stability is not certified.',
'CED-23':'Diagnostic alternative scaffold: TE/CDK1 missing with solubility120uM/half-life38min. Deposited ligand shares tautomer-normalized parent, not identical supplied state; structure is no activity prediction.',
'CED-24':'Reserve: unmeasured alternate scaffold but solubility3uM vs120uM in23; no proven shared pose or cellular activity.'}
def run(out):
 inv=readcsv(out/'analysis/inventory_identity.csv');menu={x['experiment_id']:x for x in readcsv(ROOT/'inputs/data/experiment_menu.csv')};ro={r['compound_id']:r for r in readcsv(ROOT/'inputs/data/program_readouts.csv')};selections=[];chosen={}
 for r in inv:
  cid=r['compound_id'];s=r['identity_status'];status=s if s in ('duplicate','quarantined') else ('advance' if cid in ADVANCE else 'diagnostic' if cid in DIAGNOSTIC else 'reserve')
  selections.append({'request_id':r['request_id'],'compound_id':cid,'status':status,'reason':r['identity_reason'] if s!='valid' else REASONS[cid]})
  if status in ('advance','diagnostic'):chosen[cid]=r
 experiments=[]
 for cid in ADVANCE+DIAGNOSTIC:
  for eid in ['CDK2-E1','CDK1-B','ENGAGE']:
   purpose={'CDK2-E1':'Define CDK2/cyclinE1 biochemical potency at1mM ATP; do not substitute historical multi-kinase label.','CDK1-B':'Pair with CDK2-E1 at1mM ATP to measure assay-specific biochemical CDK1/CDK2 window.','ENGAGE':'Measure or confirm cellular CDK2 engagement; use alongside matched biochemistry, not a cross-assay selectivity ratio.'}[eid]
   if cid=='CED-08' and eid=='ENGAGE':purpose+=' Two separately prepared concentration-response experiments in parallel; discordance remains unresolved rather than averaged into a pass.'
   if cid=='CED-03':purpose+=' Start only after VIAL-03 passes IDENTITY.'
   experiments.append({'compound_id':cid,'experiment_id':eid,'replicates':2 if (cid,eid)==('CED-08','ENGAGE') else 1,'purpose':purpose})
 experiments.append({'compound_id':'CED-03','experiment_id':'IDENTITY','replicates':1,'purpose':'Confirm VIAL-03 identity and purity against its intended graph before biology; keep VIAL-CONFLICT quarantined. Does not repair attribution of historic measurements.'})
 accounting=[]
 for cid,r in chosen.items():
  es=[e for e in experiments if e['compound_id']==cid];mass=sum(float(menu[e['experiment_id']]['mass_mg'])*e['replicates'] for e in es);assaycost=sum(float(menu[e['experiment_id']]['cost_usd'])*e['replicates'] for e in es);identity_days=1 if any(e['experiment_id']=='IDENTITY' for e in es) else 0;duration=max(int(menu[e['experiment_id']]['turnaround_days']) for e in es if e['experiment_id']!='IDENTITY');arrival=int(r['delivery_days']);complete=arrival+identity_days+duration
  assert float(r['purity_pct'])>=90 or identity_days==1
  assert mass<=float(r['mass_mg']) and complete<=8
  accounting.append({'compound_id':cid,'request_id':r['request_id'],'vial_id':r['vial_id'],'status':'advance' if cid in ADVANCE else 'diagnostic','handling_usd':float(r['handling_usd']),'assays_usd':assaycost,'total_usd':assaycost+float(r['handling_usd']),'required_mass_mg':mass,'available_mass_mg':float(r['mass_mg']),'remaining_mass_mg':float(r['mass_mg'])-mass,'arrival_day':arrival,'identity_days':identity_days,'parallel_assay_days':duration,'completion_day':complete,'cell_TE_nM':ro[cid]['cell_target_engagement_ic50_nM'],'CDK1_nM':ro[cid]['cdk1_cyclinB_ic50_nM'],'solubility_uM':ro[cid]['kinetic_solubility_uM'],'microsome_t_half_min':ro[cid]['microsomal_t_half_min']})
 total=sum(r['total_usd'] for r in accounting);assert total<=18000 and len(chosen)<=8 and len(ADVANCE)<=6 and len(DIAGNOSTIC)<=2;assert len({r['vial_id'] for r in chosen.values()})==len(chosen)
 writecsv(out/'delivery/selection.csv',selections);writecsv(out/'delivery/experiments.csv',experiments);writecsv(out/'analysis/cycle_accounting.csv',accounting)
 summary={'total_usd':total,'remaining_budget_usd':18000-total,'handling_usd':sum(r['handling_usd'] for r in accounting),'assays_usd':sum(r['assays_usd'] for r in accounting),'unique_compounds':len(chosen),'advancement_count':len(ADVANCE),'diagnostic_count':len(DIAGNOSTIC),'unique_vials':len(chosen),'experiment_combinations':len(experiments),'charged_menu_experiments':sum(e['replicates'] for e in experiments),'latest_completion_day':max(r['completion_day'] for r in accounting),'all_material_feasible':True,'schedule_assumption':'Arrival day0 is start; IDENTITY first adds1day; other experiments including independently prepared ENGAGE replicates run in parallel. No adaptive second wave or extra synthesis.','identity_failure_action':'Cancel CED03 biological work or quarantine uninterpretable results; do not replace from conflicting vial. Quoted total is conservative planned maximum.','gates':'Team working next-cycle criterion: identity-qualified TE<150nM, matched biochemical CDK1-B/CDK2-E1 IC50 ratio>=10, microsome>20min, adequate soluble exposure. Tenfold is a project decision rule, not supplied/clinical threshold or cellular safety claim.'}
 (out/'analysis/cycle_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--out',type=pathlib.Path,default=ROOT);run(a.parse_args().out.resolve())
