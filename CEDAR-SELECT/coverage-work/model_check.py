"""Independent artifact audit: no import/run of model-work/run.py, no fitting/tuning."""
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'): os.environ[key]='1'
import csv,json,math,statistics,time,resource,hashlib
from pathlib import Path
from collections import defaultdict,Counter
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr
start=time.time(); root=Path(__file__).resolve().parents[1]; res=root/'model-work/results'; out=Path(__file__).parent
def load(p):return json.loads(p.read_text())
def table(p):return list(csv.DictReader(p.open()))
src=load(root/'inputs/public_assay_context.json'); ana=load(root/'inputs/coverage/analogue_evidence.json'); summary=load(res/'summary.json')
raw=src['train']+src['development']+src['bounds_context_only']; by=defaultdict(dict)
def graph(s,stereo=True):return Chem.MolToSmiles(Chem.MolFromSmiles(s),isomericSmiles=stereo)
def checkrow(r):
 assert r['target'] in ['CDK2','CDK1']
 if r['relation']=='=': assert 'pIC50' in r and math.isfinite(r['pIC50']) and 'bound_pIC50' not in r
 else: assert r['relation']=='<=' and 'pIC50' not in r and r['bound_pIC50']==5.0 and r['primary_numeric'] is False
for r in raw:
 checkrow(r); assert r['target'] not in by[r['compound_id']];by[r['compound_id']][r['target']]=r
assert len({r['row_id'] for r in raw})==len(raw)
for cid,rr in by.items():
 assert set(rr)=={'CDK1','CDK2'}
 for key in ['partition','leakage_group_id','series_id','study_id']:
  assert len({r[key] for r in rr.values()})==1
 assert len({graph(r['smiles']) for r in rr.values()})==1
comp=table(res/'compounds.csv'); cmap={r['compound_id']:r for r in comp};assert set(cmap)==set(by)
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
fps={c:gen.GetFingerprint(Chem.MolFromSmiles(rr['CDK2']['smiles'])) for c,rr in by.items()}
collisions=[]
for i,c in enumerate(by):
 for d in list(by)[:i]:
  a,b=by[c]['CDK2'],by[d]['CDK2']
  reason=[]
  if a['leakage_group_id']==b['leakage_group_id']:reason.append('source_group')
  if graph(a['smiles'],False)==graph(b['smiles'],False):reason.append('nonstereo_graph')
  if fps[c]==fps[d]:reason.append('fingerprint')
  if reason:
   assert cmap[c]['merged_group']==cmap[d]['merged_group'];collisions.append(dict(compounds=[c,d],reasons=reason))
groups=defaultdict(set)
for c in comp:groups[c['merged_group']].add(c['partition'])
assert all(len(v)==1 for v in groups.values())
cross=table(res/'analogue_crosswalk.csv');assert len(cross)==16
for c in ana['compounds']:
 matches=[cid for cid,rr in by.items() if graph(c['smiles'])==graph(rr['CDK2']['smiles'])];assert len(matches)==1
 for o in ana['observations']:
  if o['compound_id']!=c['compound_id']:continue
  r=by[matches[0]][o['target']];assert o['original_unit']=='uM' and o['original_relation']=='='
  assert abs(6-math.log10(float(o['original_value']))-r['pIC50'])<1e-10
  assert any(x['context_row_id']==r['row_id'] and x['analogue_observation_id']==o['observation_id'] for x in cross)
def label(cid,endpoint):
 rr=by[cid]
 if endpoint=='log_selectivity':return rr['CDK2']['pIC50']-rr['CDK1']['pIC50']
 return rr[endpoint]['pIC50']
folds=table(res/'training_folds.csv');preds=table(res/'development_predictions.csv');metrics={};fold_checks={}
for ep,s in summary['endpoint_results'].items():
 eligible={c for c in by if ep=='CDK2' or by[c]['CDK1']['relation']=='='}
 train={c for c in eligible if by[c]['CDK2']['partition']=='train'};dev=eligible-train
 f=[r for r in folds if r['endpoint']==ep];assert len(f)==len(train) and {r['compound_id'] for r in f}==train
 foldby={r['compound_id']:r['fold'] for r in f}
 for c in train:
  for d in train:
   if cmap[c]['merged_group']==cmap[d]['merged_group']:assert foldby[c]==foldby[d]
 assert set(cmap[c]['merged_group'] for c in train).isdisjoint(cmap[c]['merged_group'] for c in dev)
 fold_checks[ep]=dict(train=len(train),development=len(dev),fold_counts=dict(Counter(r['fold'] for r in f)))
 metrics[ep]={}
 for model in ['mean','3nn','ridge']:
  pp=[r for r in preds if r['endpoint']==ep and r['model']==model];assert len(pp)==len(dev) and {r['compound_id'] for r in pp}==dev
  y=[];p=[]
  for r in pp:
   cid=r['compound_id'];yy=label(cid,ep);v=float(r['prediction']);assert abs(float(r['observed'])-yy)<1e-12
   assert abs(float(r['error'])-(v-yy))<1e-12
   assert (r['selected_on_training']=='True')==(model==s['selected'])
   sims={c:DataStructs.TanimotoSimilarity(fps[cid],fps[c]) for c in train}
   assert abs(max(sims.values())-float(r['nearest_training_tanimoto']))<1e-12
   if model=='mean':assert abs(v-statistics.mean(label(c,ep) for c in train))<1e-12
   if model=='3nn':
    # Stable source ordering independently retains equal-similarity neighbor choice.
    near=sorted([c for c in by if c in train],key=lambda c:-sims[c])[:3]
    assert abs(v-statistics.mean(label(c,ep) for c in near))<1e-12
   y.append(yy);p.append(v)
  m=dict(n=len(y),MAE=statistics.mean(abs(a-b) for a,b in zip(y,p)),RMSE=math.sqrt(statistics.mean((a-b)**2 for a,b in zip(y,p))),bias_prediction_minus_observed=statistics.mean(b-a for a,b in zip(y,p)),spearman=float(spearmanr(y,p).statistic) if max(p)-min(p)>1e-12 else None)
  for k,v in m.items():
   if v is None:assert s['development'][model][k] is None
   else:assert abs(v-s['development'][model][k])<1e-10
  metrics[ep][model]=m
 chosen=min(['mean','3nn','ridge'],key=lambda m:s['training_cv'][m]['MAE']);assert chosen==s['selected']
 gate=len(dev)>=8 and metrics[ep][chosen]['MAE']<=.8*metrics[ep]['mean']['MAE'] and (metrics[ep][chosen]['spearman'] or -1)>=.5
 assert gate==s['provisional_ranking_gate']
bounds=table(res/'bounds_checks.csv');assert {r['compound_id'] for r in bounds}=={r['compound_id'] for r in src['bounds_context_only']}
for r in bounds:assert float(r['bound_violation_log'])==max(0.,float(r['prediction_exact_only_model'])-5.)
for c in comp:
 rr=by[c['compound_id']];a=rr['CDK2']['pIC50'];b=rr['CDK1'].get('pIC50',rr['CDK1'].get('bound_pIC50'))
 assert abs(float(c['CDK2_nM'])-10**(9-a))<1e-7
 assert abs(float(c['ratio_or_lower_bound'])-10**(a-b))<1e-7
 assert c['ratio_relation']==('=' if rr['CDK1']['relation']=='=' else '>=')
neg={}
for name,mut in [('bound_wrong_inequality',dict(src['bounds_context_only'][0],relation='>=')),('bound_as_exact_label',dict(src['bounds_context_only'][0],pIC50=5.0))]:
 try:checkrow(mut);neg[name]=False
 except AssertionError:neg[name]=True
neg['micromolar_as_nanomolar_rejected']=abs(9-math.log10(.007)-8.154901959985743)>1
assert all(neg.values())
cluster=table(res/'cluster_predictions.csv'); cluster_counts={ep:len({r['compound_id'] for r in cluster if r['endpoint']==ep}) for ep in metrics}
report=dict(status='passed supplied-data checks with material reporting limitations',source_counts=dict(rows=len(raw),compounds=len(by),exact=sum(r['relation']=='=' for r in raw),bounds=len(bounds),merged_groups=len(groups)),partition_and_folds=fold_checks,collision_pairs=collisions,recalculated_development_metrics=metrics,negative_controls=neg,cluster_scored_compound_counts=cluster_counts,cluster_sizes=summary['structural_cluster_sizes'],limitations=['Training OOF predictions were not exported; selected model agrees with stored CV metrics but independent full CV fit was outside this audit.','CDK1 and ratio structural cluster scores cover only one compound each: dominant cluster excluded when fewer than three eligible training compounds remain. CDK2 covers all60 but dominant-cluster holdout trains on just3 compounds. These metrics do not establish scaffold extrapolation.','Bound code hardcodes <=5 and exact numeric eligibility uses presence of pIC50; correct supplied inputs checked but changed schemas require explicit assertions.','New-design core domain rule is stated in protocol but not exercised by this evaluation script; cannot claim candidate-domain acceptance from it alone.'],seconds=time.time()-start,cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,maxrss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,input_hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'model-work/run.py',root/'model-work/PROTOCOL.md',res/'summary.json',res/'development_predictions.csv']},jobs_stopped=True)
(out/'model_validation.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({k:report[k] for k in ['status','source_counts','cluster_scored_compound_counts','cluster_sizes','seconds','maxrss_KiB']},indent=2))
