"""Replay supplied coverage records; no external services or new measurements."""
import csv, hashlib, json, math, os, resource, time
from pathlib import Path
from rdkit import Chem, rdBase
t0=time.monotonic()
root=Path(__file__).resolve().parents[1]; out=Path(__file__).resolve().parent
def read(p): return json.loads((root/p).read_text())
def write(name,obj): (out/name).write_text(json.dumps(obj,indent=2)+'\n')
a=read('inputs/coverage/analogue_evidence.json'); ctx=read('inputs/public_assay_context.json')
b=read('inputs/coverage/entity_bindings.json'); w=read('inputs/coverage/work_orders.json'); r=read('inputs/coverage/retrieval_archive.json')
def canonical(s):
 m=Chem.MolFromSmiles(s); assert m is not None
 return Chem.MolToSmiles(m,isomericSmiles=True)
context=[x for k,v in ctx.items() if isinstance(v,list) for x in v if isinstance(x,dict) and 'smiles' in x]
rows=[]; lookup={}
for c in a['compounds']:
 obs={x['target']:x for x in a['observations'] if x['compound_id']==c['compound_id']}
 assert len(obs)==2
 for x in obs.values():
  assert x['original_unit']=='uM' and x['original_relation']=='='
  assert abs((6-math.log10(float(x['original_value'])))-x['pIC50'])<1e-10
 row=dict(c,canonical_isomeric_smiles=canonical(c['smiles']))
 row['observations']=obs
 row['CDK2_nM']=float(obs['CDK2']['original_value'])*1000
 row['CDK1_nM']=float(obs['CDK1']['original_value'])*1000
 row['ratio_CDK1_over_CDK2']=row['CDK1_nM']/row['CDK2_nM']
 row['matched_context_rows']=[x for x in context if canonical(x['smiles'])==canonical(c['smiles'])]
 row['overlap_label_checks']=[dict(row_id=x['row_id'],target=x['target'],absolute_pIC50_difference=abs(x['pIC50']-obs[x['target']]['pIC50'])) for x in row['matched_context_rows']]
 assert all(x['absolute_pIC50_difference']<1e-10 for x in row['overlap_label_checks'])
 rows.append(row); lookup[c['public_source_identifiers']['example']]=row
write('analogue_reconciled.json',rows)
pairs=[]
for before,after,edit in [(12,23,'methoxy to ethoxy'),(23,39,'terminal ethoxy H to F'),(39,41,'terminal fluoroethoxy H to F'),(41,56,'terminal difluoroethoxy H to F'),(77,12,'hydroxy O-methylation'),(69,70,'side-chain stereochemical inversion')]:
 x,y=lookup[before],lookup[after]
 pairs.append(dict(before_example=before,after_example=after,edit=edit,before_id=x['compound_id'],after_id=y['compound_id'],delta_pIC50_CDK2=y['observations']['CDK2']['pIC50']-x['observations']['CDK2']['pIC50'],delta_pIC50_CDK1=y['observations']['CDK1']['pIC50']-x['observations']['CDK1']['pIC50'],selectivity_ratio_before=x['ratio_CDK1_over_CDK2'],selectivity_ratio_after=y['ratio_CDK1_over_CDK2'],observation_ids=[v['observation_id'] for z in [x,y] for v in z['observations'].values()]))
write('local_pairs.json',pairs)
for x in b['bindings']:
 assert canonical(x['smiles'])==canonical(next(c['smiles'] for c in a['compounds'] if c['compound_id']==x['compound_id']))
ob={}
for f in w['financial_records']:
 if f['obligation_id'] in ob: assert (f['state'],f['amount'])==(ob[f['obligation_id']]['state'],ob[f['obligation_id']]['amount'])
 ob[f['obligation_id']]=f
inv={}
for v in w['inventory']:
 if v['vial_id'] in inv: assert all(v[k]==inv[v['vial_id']][k] for k in ['compound_slot','usable_mg','reserved_mg','available_business_day'])
 inv[v['vial_id']]=v
balance=dict(earmark=10000,parent_tranche=120000,paid=sum(x['amount'] for x in ob.values() if x['state']=='paid'),reserved=sum(x['amount'] for x in ob.values() if x['state']=='reserved'),obligations=list(ob.values()),inventory=list(inv.values()))
balance.update(uncommitted_before=10000-balance['paid']-balance['reserved'],proposed=['O1','O2'],proposed_incremental=3400,uncommitted_after=2600,total_earmark_committed_after=7400,inventory_free_after={'V1':4,'V2':4,'V3':2},schedule=[{'work':'W2/CDK1','lane':1,'start':0,'finish':3,'already_authorized':True},{'work':'O2','lane':2,'start':0,'finish':3},{'work':'O1','lane':2,'start':3,'finish':5}])
assert balance['uncommitted_before']==6000
write('budget_material_schedule.json',balance)
classes=[('R1','qualified index only','Intent/execution/target/endpoint agree; primary measurement absent.'),('R2','unsupported executed filter and inapplicable returned target','target is not target_accession; returned P06493 fails intended P24941.'),('R3','wrong typed entity','catalog:assay:771 is not catalog:compound:771; reissue compound query.'),('R4','partial','Only page 1 of 2 supplied; cannot claim complete coverage.'),('R5','valid empty','Completed supported compound 999 lookup; no universal absence/inactivity inference.'),('R6','failed','HTTP 503 is not an empty biological result.'),('R7','valid empty child lookup','Documented parent relationship permits bounded parent mechanism lookup R8.'),('R8','parent annotation only','Parent-indexed annotation cannot transfer numerical potency or establish actual salt identity.'),('R9','inapplicable endpoint','percent inhibition at 10 uM does not satisfy IC50 despite matching executed endpoint filter.')]
write('retrieval_audit.json',[dict(response_id=i,status=s,decision=d,original=next(x for x in r['responses'] if x['id']==i)) for i,s,d in classes])
write('resource_log.json',dict(start_utc_epoch=time.time()-(time.monotonic()-t0),wall_seconds=time.monotonic()-t0,cpu_seconds=resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,maxrss_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,rdkit=rdBase.rdkitVersion,threads={k:os.environ.get(k) for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']},input_sha256={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'inputs/public_assay_context.json',*sorted((root/'inputs/coverage').glob('*.json'))]},completed=True))
print(json.dumps([{'example':x['public_source_identifiers']['example'],'CDK2_nM':x['CDK2_nM'],'CDK1_nM':x['CDK1_nM'],'ratio':x['ratio_CDK1_over_CDK2'],'overlap_rows':len(x['matched_context_rows'])} for x in rows],indent=2))
