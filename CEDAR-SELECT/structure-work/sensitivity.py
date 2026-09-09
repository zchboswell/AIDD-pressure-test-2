import json,time,resource
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SVDSuperimposer import SVDSuperimposer
from pathlib import Path
start=time.time();root=Path('CEDAR-SELECT/structure-work');r=json.loads((root/'analysis.json').read_text());result={'comparisons':[],'sulfonamide_Asp86':[]}
def atoms(acc):
 d=MMCIF2Dict('CEDAR-SELECT/inputs/structural_sources/'+acc+'.cif');ks=[k for k in d if k.startswith('_atom_site.')]; aa=[dict(zip([k.split('.')[1] for k in ks],v)) for v in zip(*(d[k] for k in ks))]
 for a in aa:a['xyz']=np.array([float(a['Cartn_'+c]) for c in 'xyz'])
 return aa
allats={a:atoms(a) for a in ('5NEV','5LQF')}
for c in r['comparisons'][2:]:
 ac,ch=c['first'].split(':');bc,dh=c['second'].split(':');ats=allats[ac];bts=allats[bc]
 a={int(x['auth_seq_id']):x for x in ats if x['auth_asym_id']==ch and x['auth_atom_id']=='CA' and x['group_PDB']=='ATOM'};b={int(x['auth_seq_id']):x for x in bts if x['auth_asym_id']==dh and x['auth_atom_id']=='CA' and x['group_PDB']=='ATOM'}
 mapping={int(i):j for i,j in c['residue_mapping'].items()};fit=[(i,j) for i,j in mapping.items() if 1<=i<=90 and i not in range(9,20) and a[i]['auth_comp_id']==b[j]['auth_comp_id']]
 xx=np.array([a[i]['xyz'] for i,j in fit]); yy=np.array([b[j]['xyz'] for i,j in fit]);sv=SVDSuperimposer();sv.set(yy,xx);sv.run();rot,tran=sv.get_rotran()
 aa=next(x for x in ats if x['auth_asym_id']==ch and x['auth_seq_id']=='15' and x['auth_atom_id']=='OH');bb=next(x for x in bts if x['auth_asym_id']==dh and x['auth_seq_id']=='15' and x['auth_atom_id']=='OH')
 sidecontacts=[]
 for x in ats:
  if x['auth_asym_id']!=ch or x['group_PDB']!='ATOM' or x['type_symbol']=='H' or x['auth_atom_id'] in ('CA','C','N','O'):continue
  lat=[l for l in ats if l['auth_asym_id']==ch and l['label_comp_id']=='72L']
  dd=min(float(np.linalg.norm(x['xyz']-l['xyz'])) for l in lat)
  if dd<=4:
   i=int(x['auth_seq_id']);j=mapping.get(i);sidecontacts.append({'CDK2_residue':i,'CDK2_name':x['auth_comp_id'],'CDK1_residue':j,'CDK1_name':None if j not in b else b[j]['auth_comp_id'],'atom':x['auth_atom_id'],'distance_A':dd})
 result['comparisons'].append({'first':c['first'],'second':c['second'],'N_lobe_fit_n':len(fit),'N_lobe_CA_RMSD_A':sv.get_rms(),'Tyr15_OH_displacement_A':float(np.linalg.norm(aa['xyz']@rot+tran-bb['xyz'])),'CDK2_sidechain_contacts':sidecontacts})
for acc,lig,chs,natom in [('5NEV','72L',['A','C'],'N6'),('5LQF','4SP',['A','D'],'N26')]:
 for ch in chs:
  aa=next(x for x in allats[acc] if x['auth_asym_id']==ch and x['label_comp_id']==lig and x['auth_atom_id']==natom)
  for at in ('OD1','OD2','N'):
   bb=next(x for x in allats[acc] if x['auth_asym_id']==ch and x['auth_seq_id']=='86' and x['auth_atom_id']==at)
   result['sulfonamide_Asp86'].append({'structure':acc,'chain':ch,'ligand_N':natom,'Asp86_atom':at,'distance_A':float(np.linalg.norm(aa['xyz']-bb['xyz']))})
result['runtime']={'elapsed_seconds':time.time()-start,'cpu_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,'maxrss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'threads':1,'stop':'normal exit'}
(root/'sensitivity.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
