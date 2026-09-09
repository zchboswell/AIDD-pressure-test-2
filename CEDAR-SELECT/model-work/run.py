"""Reproducible local evaluation; scientific contract in PROTOCOL.md."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']: os.environ[k]='1'
import json,csv,time,resource,sys,hashlib
from pathlib import Path
import numpy as np
from rdkit import Chem,DataStructs,rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
from sklearn.model_selection import GroupKFold
from sklearn.kernel_ridge import KernelRidge
from scipy.stats import spearmanr
START=time.time(); ROOT=Path(__file__).resolve().parents[1]; OUT=Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'model-work/results'; OUT.mkdir(parents=True,exist_ok=True)
def save(name,obj): (OUT/name).write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n')
def csvout(name,rows):
 if not rows:return
 with (OUT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
j=json.loads((ROOT/'inputs/public_assay_context.json').read_text()); raw=j['train']+j['development']+j['bounds_context_only']; by={};cur=[]
for r in raw:
 m=Chem.MolFromSmiles(r['smiles']); assert m is not None
 s=Chem.MolToSmiles(m,True); cid=r['compound_id']; d=by.setdefault(cid,{'id':cid,'smiles':s,'partition':r['partition'],'group':r['leakage_group_id']})
 assert d['smiles']==s and d['partition']==r['partition'] and r['target'] not in d
 d[r['target']]=r
 cur.append(dict(row_id=r['row_id'],compound_id=cid,target=r['target'],partition=r['partition'],source_relation=r['relation'],pIC50=r.get('pIC50',''),bound_pIC50=r.get('bound_pIC50',''),status='exact_regression' if r['relation']=='=' else 'bound_context',canonical_isomeric_smiles=s))
D=list(by.values());N=len(D); mol=[Chem.MolFromSmiles(d['smiles']) for d in D];gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True);fp=[gen.GetFingerprint(m) for m in mol];X=np.array([list(f) for f in fp],dtype=float)
inter=X@X.T; S=inter/(X.sum(1)[:,None]+X.sum(1)[None,:]-inter)
parent=list(range(N))
def find(i):
 while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
 return i
def union(i,k):parent[find(k)]=find(i)
noniso=[Chem.MolToSmiles(m,isomericSmiles=False) for m in mol]
for i in range(N):
 for k in range(i):
  if D[i]['group']==D[k]['group'] or noniso[i]==noniso[k] or np.array_equal(X[i],X[k]):union(i,k)
groups=np.array([find(i) for i in range(N)]);mixed={g for g in groups if len({D[i]['partition'] for i in range(N) if groups[i]==g})>1};tr=np.array([i for i,d in enumerate(D) if d['partition']=='train' and groups[i] not in mixed]);te=np.array([i for i,d in enumerate(D) if d['partition']=='development' and groups[i] not in mixed]);assert not set(groups[tr])&set(groups[te])
# Construct exact endpoints without converting bounds into point labels.
Y={t:np.array([d[t].get('pIC50',np.nan) for d in D]) for t in ['CDK2','CDK1']};Y['log_selectivity']=Y['CDK2']-Y['CDK1']
MODELS=['mean','3nn','ridge']
def pred(model,train,test,y):
 if model=='mean':return np.repeat(np.mean(y[train]),len(test))
 if model=='3nn':return np.array([np.mean(y[train[np.argsort(-S[i,train],kind='stable')[:3]]]) for i in test])
 mu=np.mean(y[train]);alpha=np.linalg.solve(S[np.ix_(train,train)]+np.eye(len(train)),y[train]-mu);return mu+S[np.ix_(test,train)]@alpha

def metric(y,p):
 rho=float(spearmanr(y,p).statistic) if len(y)>2 and np.std(p)>1e-12 else None
 return {'n':len(y),'MAE':float(np.mean(abs(y-p))),'RMSE':float(np.sqrt(np.mean((y-p)**2))),'bias_prediction_minus_observed':float(np.mean(p-y)),'spearman':rho}
summary={};predrows=[];foldrows=[];oofall={};selected={};bands={}
for endpoint,y in Y.items():
 a=tr[np.isfinite(y[tr])]; b=te[np.isfinite(y[te])];cv=GroupKFold(n_splits=min(5,len(set(groups[a]))));oof={m:np.full(N,np.nan) for m in MODELS}
 for fold,(ai,bi) in enumerate(cv.split(a,groups=groups[a])):
  train,test=a[ai],a[bi];assert not set(groups[train])&set(groups[test])
  for i in test:foldrows.append({'endpoint':endpoint,'fold':fold,'compound_id':D[i]['id'],'group':int(groups[i])})
  for model in MODELS:oof[model][test]=pred(model,train,test,y)
 cvmetrics={m:metric(y[a],oof[m][a]) for m in MODELS};chosen=min(MODELS,key=lambda m:cvmetrics[m]['MAE']);selected[endpoint]=chosen
 band=float(np.quantile(abs(y[a]-oof[chosen][a]),.9,method='higher'));bands[endpoint]=band;dev={}
 for model in MODELS:
  p=pred(model,a,b,y);dev[model]=metric(y[b],p)
  for i,v in zip(b,p):predrows.append({'endpoint':endpoint,'model':model,'selected_on_training':model==chosen,'compound_id':D[i]['id'],'observed':float(y[i]),'prediction':float(v),'error':float(v-y[i]),'nearest_training_tanimoto':float(max(S[i,a])),'descriptive_error_band':band if model==chosen else ''})
 accept=bool(len(b)>=8 and dev[chosen]['MAE']<=.8*dev['mean']['MAE'] and (dev[chosen]['spearman'] or -1)>=.5)
 summary[endpoint]={'train_n':len(a),'development_n':len(b),'training_cv':cvmetrics,'selected':chosen,'development':dev,'provisional_ranking_gate':accept,'descriptive_90pct_abs_error':band};oofall[endpoint]=oof
# Structural cluster stress test; no label-informed clustering.
dists=[float(1-S[i,k]) for i in range(N) for k in range(i)];clusters=Butina.ClusterData(dists,N,.35,isDistData=True,reordering=True)
clusterparent=list(range(N))
def cf(i):
 while clusterparent[i]!=i:i=clusterparent[i]
 return i
def cu(i,k):clusterparent[cf(k)]=cf(i)
for cl in clusters:
 for i in cl[1:]:cu(cl[0],i)
for i in range(N):
 for k in range(i):
  if groups[i]==groups[k]:cu(i,k)
cg=np.array([cf(i) for i in range(N)]); stress={};stressrows=[]
if len(set(cg))>=3:
 for endpoint,y in Y.items():
  ids=np.where(np.isfinite(y))[0];oof={m:np.full(N,np.nan) for m in MODELS}
  for g in sorted(set(cg[ids])):
   a=ids[cg[ids]!=g];b=ids[cg[ids]==g]
   if len(a)<3:continue
   for model in MODELS:
    p=pred(model,a,b,y);oof[model][b]=p
    for i,v in zip(b,p):stressrows.append({'endpoint':endpoint,'model':model,'compound_id':D[i]['id'],'cluster':int(g),'observed':float(y[i]),'prediction':float(v),'nearest_train_similarity':float(max(S[i,a]))})
  stress[endpoint]={m:metric(y[np.isfinite(oof[m])],oof[m][np.isfinite(oof[m])]) for m in MODELS}
# Group-preserving permutation null for fixed mean versus fixed 3NN train CV.
y=Y['CDK2'];a=tr; cv=list(GroupKFold(n_splits=5).split(a,groups=groups[a]));rng=np.random.default_rng(20260909);gr={g:np.array([i for i in a if groups[i]==g]) for g in set(groups[a])};strata={n:[g for g in gr if len(gr[g])==n] for n in set(map(len,gr.values()))};null=[]
for rep in range(100):
 yp=y.copy()
 for n,gg in strata.items():
  shuffled=rng.permutation(gg)
  for dest,src in zip(gg,shuffled):yp[gr[dest]]=y[gr[src]]
 errs={m:[] for m in ['mean','3nn']}
 for ai,bi in cv:
  aa,bb=a[ai],a[bi]
  for m in errs:errs[m].extend(abs(yp[bb]-pred(m,aa,bb,yp)))
 null.append(float(np.mean(errs['mean'])-np.mean(errs['3nn'])))
observed=summary['CDK2']['training_cv']['mean']['MAE']-summary['CDK2']['training_cv']['3nn']['MAE']
# Bounds retained as inequalities and scored only for violations, not point errors.
boundrows=[]
for i,d in enumerate(D):
 if d['CDK1']['relation']=='=':continue
 y=Y['CDK1'];a=tr[np.isfinite(y[tr])];p=float(pred(selected['CDK1'],a,np.array([i]),y)[0]);boundrows.append({'compound_id':d['id'],'partition':d['partition'],'CDK1_pIC50_upper_bound':5.,'prediction_exact_only_model':p,'bound_violation_log':max(0.,p-5.),'interpretation':'bound retained; exact-only model selection bias'})
# Crosswalk is measurement reuse, never additional rows.
ana=json.loads((ROOT/'inputs/coverage/analogue_evidence.json').read_text());cross=[]
for c in ana['compounds']:
 s=Chem.MolToSmiles(Chem.MolFromSmiles(c['smiles']),True);matches=[d for d in D if d['smiles']==s];assert len(matches)==1
 d=matches[0]
 for r in [r for r in ana['observations'] if r['compound_id']==c['compound_id']]:
  rr=d[r['target']];assert abs(r['pIC50']-rr['pIC50'])<1e-10 and r['p_relation']==rr['relation']
  cross.append({'analogue_compound_id':c['compound_id'],'example':c['public_source_identifiers']['example'],'analogue_observation_id':r['observation_id'],'context_compound_id':d['id'],'context_row_id':rr['row_id'],'target':r['target'],'action':'same graph and label; count once'})
# Independent numeric implementations and chemical invariants.
assert abs((9-np.log10(7))-8.154901959985743)<1e-12
assert abs(10**(8.154901959985743-5.920818753952375)-1200/7)<1e-9
assert 9-np.log10(10000)==5 and 9-np.log10(20000)<5
assert np.max(abs(S-np.array([[DataStructs.TanimotoSimilarity(f,g) for g in fp] for f in fp])))<1e-12
assert np.max(abs(S-S.T))<1e-12 and np.allclose(np.diag(S),1)
y=Y['CDK2'];a=tr;b=te;mu=np.mean(y[a]);kp=KernelRidge(alpha=1,kernel='precomputed').fit(S[np.ix_(a,a)],y[a]-mu).predict(S[np.ix_(b,a)])+mu;assert np.max(abs(kp-pred('ridge',a,b,y)))<1e-10
compoundrows=[]
for i,d in enumerate(D):
 r1=d['CDK1'];v1=r1.get('pIC50',r1.get('bound_pIC50'));v2=d['CDK2']['pIC50'];compoundrows.append({'compound_id':d['id'],'partition':d['partition'],'leakage_group':d['group'],'merged_group':int(groups[i]),'structural_cluster':int(cg[i]),'smiles':d['smiles'],'CDK2_nM':10**(9-v2),'CDK1_nM_or_lower_bound':10**(9-v1),'CDK1_IC50_relation':'=' if r1['relation']=='=' else '>=','ratio_or_lower_bound':10**(v2-v1),'ratio_relation':'=' if r1['relation']=='=' else '>=','joint_biochemical_gate':v2>=7 and v2-v1>=1})
save('summary.json',{'endpoint_results':summary,'joint_model_led_selection':summary['CDK2']['provisional_ranking_gate'] and summary['log_selectivity']['provisional_ranking_gate'],'source_rows':len(raw),'compounds':N,'merged_identity_groups':len(set(groups)),'mixed_partition_groups_excluded':list(map(int,mixed)),'structural_cluster_sizes':sorted([int(sum(cg==g)) for g in set(cg)],reverse=True),'cluster_stress':stress,'permutation_CDK2':{'repeats':100,'observed_MAE_improvement_3nn_vs_mean':observed,'null_improvements':null,'one_sided_empirical_p':(1+sum(v>=observed for v in null))/101},'independent_checks':'passed analytic transforms, bounds, RDKit similarity, sklearn ridge, group boundaries','scope':'nonblind same-public-study development; no unseen-study, cellular or full A1 claims'})
for name,rs in [('curation.csv',cur),('compounds.csv',compoundrows),('development_predictions.csv',predrows),('training_folds.csv',foldrows),('cluster_predictions.csv',stressrows),('bounds_checks.csv',boundrows),('analogue_crosswalk.csv',cross)]:csvout(name,rs)
save('runtime.json',{'seconds':time.time()-START,'maxrss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'python':sys.version,'rdkit':rdBase.rdkitVersion,'numpy':np.__version__,'threads':1,'input_sha256':hashlib.sha256((ROOT/'inputs/public_assay_context.json').read_bytes()).hexdigest(),'code_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()})
print(json.dumps({k:{'selected':v['selected'],'development':v['development'],'accepted':v['provisional_ranking_gate']} for k,v in summary.items()},indent=2));print('groups',len(set(groups)),'clusters',sorted([int(sum(cg==g)) for g in set(cg)],reverse=True),'runtime',time.time()-START)
