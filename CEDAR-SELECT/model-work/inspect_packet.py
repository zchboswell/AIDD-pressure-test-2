import json,collections
j=json.load(open('CEDAR-SELECT/inputs/public_assay_context.json'))
rows=j['train']+j['development']+j['bounds_context_only'];d={}
for r in rows:d.setdefault(r['compound_id'],{})[r['target']]=r
for p in ['train','development']:
 ds=[v for v in d.values() if v['CDK2']['partition']==p]
 print(p,'compounds',len(ds),'paired exact',sum(v['CDK1']['relation']=='=' for v in ds))
 for t in ['CDK2','CDK1']:
  vs=[v[t]['pIC50'] for v in ds if v[t]['relation']=='='];print(t,min(vs),max(vs))
 print('biochem_gate',sum(v['CDK2']['pIC50']>=7 and v['CDK2']['pIC50']-v['CDK1'].get('pIC50',v['CDK1'].get('bound_pIC50'))>=1 for v in ds))
print('unique labels',len(rows),'groups',len(set(r['leakage_group_id'] for r in rows)))
