"""Diagnostics only: failing model gate prevents numerical advancement claims."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
import json,csv,time,resource
from pathlib import Path
import numpy as np
from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
ROOT=Path(__file__).resolve().parents[1];start=time.time()
j=json.loads((ROOT/'inputs/public_assay_context.json').read_text());rows=j['train'];D={}
for r in rows:D.setdefault(r['compound_id'],{'smiles':r['smiles'],'id':r['compound_id']})[r['target']]=r['pIC50']
ds=list(D.values());gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True);fp=[gen.GetFingerprint(Chem.MolFromSmiles(d['smiles'])) for d in ds];S=np.array([[DataStructs.TanimotoSimilarity(f,g) for g in fp] for f in fp]);summary=json.loads((ROOT/'model-work/results/summary.json').read_text());designs=json.loads((ROOT/'design-work/designs.json').read_text());out=[]
core=Chem.MolFromSmiles('O=C1N([C@@H]2CCC[C@@H](O)C2)c2nc(Nc3c[nH]nc3)ncc2C12CC2')
assert core is not None and not Chem.MolFromSmiles('c1ccccc1').HasSubstructMatch(core,useChirality=True)
for d in designs:
 f=gen.GetFingerprint(Chem.MolFromSmiles(d['smiles']));sim=np.array(DataStructs.BulkTanimotoSimilarity(f,fp));r={'design_id':d['design_id'],'smiles':d['canonical_smiles'],'status':'historical_observation_not_prediction' if d['supplied_matches'] else 'abstain_from_model_led_selection; numerical_diagnostics_only','nearest_training_id':ds[int(sim.argmax())]['id'],'nearest_training_tanimoto':float(max(sim)),'similarity_domain':bool(max(sim)>=.65),'same_stereochemical_core':Chem.MolFromSmiles(d['smiles']).HasSubstructMatch(core,useChirality=True),'uncertainty':'unavailable calibrated uncertainty; model gate failed'}
 for target in ['CDK2','CDK1','log_selectivity']:
  y=np.array([x.get(target,np.nan) if target!='log_selectivity' else x.get('CDK2',np.nan)-x.get('CDK1',np.nan) for x in ds]);ids=np.where(np.isfinite(y))[0];mu=y[ids].mean();chosen=summary['endpoint_results'][target]['selected'];nn=y[ids[np.argsort(-sim[ids],kind='stable')[:3]]].mean();ridge=mu+sim[ids]@np.linalg.solve(S[np.ix_(ids,ids)]+np.eye(len(ids)),y[ids]-mu);p={'mean':mu,'3nn':nn,'ridge':ridge}
  r[target+'_diagnostic']=float(p[chosen]);r[target+'_model_range_low']=float(min(p.values()));r[target+'_model_range_high']=float(max(p.values()));r[target+'_training_error_band']=summary['endpoint_results'][target]['descriptive_90pct_abs_error']
 assert r['same_stereochemical_core']
 out.append(r)
(ROOT/'model-work/results/design_diagnostics.json').write_text(json.dumps(out,indent=2)+'\n')
(ROOT/'model-work/results/design_runtime.json').write_text(json.dumps({'seconds':time.time()-start,'maxrss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'threads':1},indent=2)+'\n')
print(json.dumps(out,indent=2))
