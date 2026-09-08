"""Frozen small audit; no deployment fit or cellular potency extrapolation."""
import argparse,collections,csv,hashlib,json,math,pathlib,subprocess,sys
import numpy as np
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator,Descriptors
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold
from sklearn.metrics import mean_absolute_error,mean_squared_error,r2_score
from curate import readcsv,writecsv
ROOT=pathlib.Path(__file__).resolve().parents[1]
def metrics(y,p,w):return {'mae':float(mean_absolute_error(y,p,sample_weight=w)),'rmse':float(mean_squared_error(y,p,sample_weight=w)**.5),'r2':float(r2_score(y,p,sample_weight=w))}
def run(out):
 out.joinpath('analysis').mkdir(parents=True,exist_ok=True)
 legacy=subprocess.run([sys.executable,str(ROOT/'inputs/legacy/vendor_model.py'),str(ROOT/'inputs')],capture_output=True,text=True)
 (out/'analysis/vendor_original_stdout.txt').write_text(legacy.stdout);(out/'analysis/vendor_original_stderr.txt').write_text(legacy.stderr)
 rows=[r for r in readcsv(out/'analysis/curation_details.csv') if r['status']=='eligible'];features={r['row_id']:r for r in readcsv(ROOT/'inputs/data/vendor_features.csv')}
 y=np.array([float(r['p_activity']) for r in rows]);chem=[r['chemical_group'] for r in rows];counts=collections.Counter(chem);w=np.array([1/counts[c] for c in chem]);mols=[Chem.MolFromSmiles(r['canonical_isomeric_smiles']) for r in rows]
 gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True);fps=[gen.GetFingerprint(m) for m in mols]
 desc=np.array([[float(features[r['row_id']][k]) for k in ['molecular_weight','logp','heavy_atoms']] for r in rows]);leaky=np.column_stack([desc,[float(features[r['row_id']]['vendor_activity_index']) for r in rows]])
 recomputed=np.array([[Descriptors.MolWt(m),Descriptors.MolLogP(m),m.GetNumHeavyAtoms()] for m in mols]);descriptor_max_abs=np.max(abs(recomputed-desc),axis=0).tolist()
 # Union scaffold, chemical, and chosen representation collisions before splitting.
 parent=list(range(len(rows)))
 def find(i):
  while i!=parent[i]:parent[i]=parent[parent[i]];i=parent[i]
  return i
 def union(i,j):parent[find(j)]=find(i)
 for keys in ([r['murcko_scaffold'] for r in rows],chem,[fp.ToBitString() for fp in fps]):
  prior={}
  for i,k in enumerate(keys):
   if k in prior:union(i,prior[k])
   else:prior[k]=i
 groups=np.array([find(i) for i in range(len(rows))]);n_groups=len(set(groups));assert n_groups>=2
 splits=list(GroupKFold(n_splits=min(5,n_groups)).split(desc,y,groups));oof={k:np.empty(len(rows)) for k in ['mean','nearest_neighbor','rf_descriptors','rf_retrospective_index']};foldrows=[];partition=[]
 for fold,(train,test) in enumerate(splits):
  assert not set(groups[train])&set(groups[test]);assert not set(np.array(chem)[train])&set(np.array(chem)[test])
  assert not set(fps[i].ToBitString() for i in train)&set(fps[i].ToBitString() for i in test)
  oof['mean'][test]=np.average(y[train],weights=w[train])
  trainchem=collections.defaultdict(list)
  for i in train:trainchem[chem[i]].append(i)
  reps=[v[0] for v in trainchem.values()]
  for i in test:
   similarities=DataStructs.BulkTanimotoSimilarity(fps[i],[fps[j] for j in reps]);best=int(np.argmax(similarities));j=reps[best];oof['nearest_neighbor'][i]=np.mean(y[trainchem[chem[j]]])
   partition.append({'row_id':rows[i]['row_id'],'chemical_group':chem[i],'scaffold_group':str(groups[i]),'fold':fold,'max_train_tanimoto':max(similarities),'nearest_train_row':rows[j]['row_id']})
  for key,x in [('rf_descriptors',desc),('rf_retrospective_index',leaky)]:
   model=RandomForestRegressor(n_estimators=80,random_state=7,n_jobs=1);model.fit(x[train],y[train],sample_weight=w[train]);oof[key][test]=model.predict(x[test])
  for key,p in oof.items():foldrows.append({'fold':fold,'method':key,'test_rows':len(test),'test_chemicals':len(set(chem[i] for i in test)),'test_scaffold_groups':len(set(groups[test])),**metrics(y[test],p[test],w[test])})
 scores={k:metrics(y,p,w) for k,p in oof.items()};predrows=[{'row_id':r['row_id'],'chemical_group':chem[i],'observed_pIC50':y[i],'chemical_weight':w[i],**{k:float(p[i]) for k,p in oof.items()}} for i,r in enumerate(rows)]
 repeats=[]
 for c,n in counts.items():
  if n>1:
   idx=[i for i,v in enumerate(chem) if v==c];repeats.append({'chemical_group':c,'row_ids':';'.join(rows[i]['row_id'] for i in idx),'n_source_records':n,'min_p':min(y[idx]),'max_p':max(y[idx]),'spread_log':float(np.ptp(y[idx]))})
 writecsv(out/'analysis/model_oof.csv',predrows);writecsv(out/'analysis/model_partitions.csv',partition);writecsv(out/'analysis/model_fold_metrics.csv',foldrows);writecsv(out/'analysis/repeated_chemistry.csv',repeats)
 inventory=readcsv(out/'analysis/inventory_identity.csv');domain=[]
 for r in inventory:
  if r['identity_status']!='valid':continue
  fp=gen.GetFingerprint(Chem.MolFromSmiles(r['canonical_isomeric_smiles']));sims=DataStructs.BulkTanimotoSimilarity(fp,fps);best=int(np.argmax(sims));matches=[rows[i] for i,c in enumerate(chem) if c==r['chemical_group']]
  domain.append({'request_id':r['request_id'],'compound_id':r['compound_id'],'max_historical_tanimoto':max(sims),'nearest_historical_row':rows[best]['row_id'],'exact_chemical_training_rows':';'.join(x['row_id'] for x in matches),'eligible_historical_observation_pIC50s':';'.join(x['p_activity'] for x in matches),'inference_role':'No cellular inference; exact historical matches are observations, not prospective validations'})
 writecsv(out/'analysis/inventory_historical_context.csv',domain)
 summary={'original_vendor_exit_code':legacy.returncode,'original_vendor_stdout':legacy.stdout.strip(),'curated_rows':len(rows),'unique_chemicals':len(counts),'scaffold_collision_groups':n_groups,'folds':len(splits),'fold_test_rows':[len(b) for a,b in splits],'metrics_chemical_weighted':scores,'descriptor_recomputation_max_absolute_difference':dict(zip(['MW','logP','heavy_atoms'],descriptor_max_abs)),'retrospective_index_correlation_with_curated_target':float(np.corrcoef(leaky[:,-1],y)[0,1]),'max_repeat_spread_log':max(float(r['spread_log']) for r in repeats),'null_control_status':'Omitted: distinct-source repeated activities disagree and row-specific cyclin context unresolved; no unique experimental unit for the originally proposed shuffle. No null-test claim.','uncertainty':'No calibrated intervals; fold spread is diagnostic not experimental precision.','deployment_decision':'Reject original vendor model for prospective prioritization. Retrospective feature unavailable at decision time and target is historical biochemical export, not current cellular endpoint. Descriptor and nearest-neighbor models are audit controls only.','evaluation_scope':'Author-visible exploratory scaffold-group evaluation; no tuning, no external data, no hidden/temporal validation claim.'}
 (out/'analysis/model_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=pathlib.Path,default=ROOT);run(ap.parse_args().out.resolve())
