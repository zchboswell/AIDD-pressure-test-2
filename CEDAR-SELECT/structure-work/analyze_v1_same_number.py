"""Immutable deposited heavy-atom contact inspection; no preparation, docking or energy inference.
Prespecified: all ligand copies; contacts <=4.0 A, polar proximity <=3.5 A;
protein comparison on common same-number identical-residue CA excluding residues 9-19;
report loop shifts without interpreting different-ligand crystal states as selectivity free energy.
"""
import json,hashlib,time,resource,math,itertools,platform
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import Bio,numpy as np
from rdkit import Chem
start=time.time(); base=Path('CEDAR-SELECT/inputs/structural_sources'); out=Path('CEDAR-SELECT/structure-work')
def rows(d,prefix):
 keys=[k for k in d if k.startswith(prefix+'.')]
 return [dict(zip([k.split('.',1)[1] for k in keys],v)) for v in zip(*(d[k] for k in keys))] if keys else []
def load(acc):
 d=MMCIF2Dict(str(base/f'{acc}.cif')); ats=rows(d,'_atom_site')
 for a in ats:a['xyz']=np.array([float(a[f'Cartn_{c}']) for c in 'xyz'])
 return d,ats
def ident(a):return ':'.join(a.get(k,'?') for k in ['auth_asym_id','auth_seq_id','auth_comp_id','auth_atom_id','label_alt_id'])
def heavy(a):return a['type_symbol'] not in ('H','D')
def dist(a,b):return float(np.linalg.norm(a['xyz']-b['xyz']))
report={'protocol':__doc__,'sources':{},'structures':{},'chemical_components':{},'comparisons':[]}
for p in base.iterdir():
 if p.is_file():report['sources'][p.name]={'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'bytes':p.stat().st_size}
meta=json.loads((base/'SOURCE_METADATA.json').read_text()); assert all(report['sources'][k]['sha256']==v['sha256'] for k,v in meta['sources'].items())
structures={}
for acc,lig in [('5NEV','72L'),('5LQF','4SP')]:
 d,ats=load(acc); structures[acc]=(d,ats)
 r={'title':d['_struct.title'],'resolution_A':d['_refine.ls_d_res_high'],'entities':rows(d,'_entity'),'polymers':rows(d,'_entity_poly'),'sequence_differences':rows(d,'_struct_ref_seq_dif'),'unobserved_residues':rows(d,'_pdbx_unobs_or_zero_occ_residues'),'unobserved_atoms':rows(d,'_pdbx_unobs_or_zero_occ_atoms'),'assemblies':rows(d,'_pdbx_struct_assembly_gen'),'revision':rows(d,'_pdbx_audit_revision_history'),'ligands':[]}
 for ch,res in sorted(set((a['auth_asym_id'],a['auth_seq_id']) for a in ats if a['label_comp_id']==lig)):
  lat=[a for a in ats if a['label_comp_id']==lig and a['auth_asym_id']==ch and a['auth_seq_id']==res and heavy(a)]
  pat=[a for a in ats if a['group_PDB']=='ATOM' and heavy(a)]
  pairs=[(dist(a,b),a,b) for a in lat for b in pat if dist(a,b)<=4]
  minima={}
  for dd,a,b in pairs:
   k=f"{b['auth_asym_id']}:{b['auth_seq_id']}:{b['auth_comp_id']}"
   if k not in minima or dd<minima[k]['distance_A']:minima[k]={'distance_A':round(dd,4),'ligand_atom':a['auth_atom_id'],'protein_atom':b['auth_atom_id']}
  # Independent distance implementation checks every reported minimum (math.dist on lists).
  assert all(abs(dd-math.dist(a['xyz'].tolist(),b['xyz'].tolist()))<1e-10 for dd,a,b in pairs)
  polar=[{'ligand_atom':ident(a),'protein_atom':ident(b),'distance_A':round(dd,4)} for dd,a,b in pairs if dd<=3.5 and a['type_symbol'] in ('N','O') and b['type_symbol'] in ('N','O')]
  waters=[{'water':ident(b),'ligand_atom':a['auth_atom_id'],'distance_A':round(dist(a,b),4)} for a in lat for b in ats if b['label_comp_id']=='HOH' and heavy(b) and dist(a,b)<=3.5]
  r['ligands'].append({'author_chain':ch,'label_chain':lat[0]['label_asym_id'],'residue':res,'heavy_atom_count':len(lat),'atom_ids':[a['auth_atom_id'] for a in lat],'occupancy_range':[min(float(a['occupancy']) for a in lat),max(float(a['occupancy']) for a in lat)],'B_range':[min(float(a['B_iso_or_equiv']) for a in lat),max(float(a['B_iso_or_equiv']) for a in lat)],'alt_ids':sorted(set(a['label_alt_id'] for a in lat)),'residue_contact_minima':minima,'polar_proximities':polar,'nearby_waters':waters})
 report['structures'][acc]=r
for ligand in ('72L','4SP'):
 d=MMCIF2Dict(str(base/f'{ligand}.cif')); atoms=rows(d,'_chem_comp_atom'); bonds=rows(d,'_chem_comp_bond'); ads={a['atom_id']:a for a in atoms}; heavyids={a['atom_id'] for a in atoms if a['type_symbol']!='H'}
 r={'name':d['_chem_comp.name'],'formula':d['_chem_comp.formula'],'formal_charge':d['_chem_comp.pdbx_formal_charge'],'descriptors':rows(d,'_pdbx_chem_comp_descriptor'),'heavy_atom_count':len(heavyids),'coordinate_bond_checks':[]}
 # Compare deposited bond lengths with independent dictionary IDEAL geometry, a diagnostic not z-score.
 for acc,(sd,sats) in structures.items():
  for lr in report['structures'][acc]['ligands']:
   if set(lr['atom_ids'])!=heavyids:continue
   lookup={a['auth_atom_id']:a for a in sats if a['auth_asym_id']==lr['author_chain'] and a['auth_seq_id']==lr['residue'] and a['label_comp_id']==ligand and heavy(a)}
   assert set(lookup)==heavyids
   vals=[]
   for b in bonds:
    a1,a2=b['atom_id_1'],b['atom_id_2']
    if a1 not in lookup or a2 not in lookup:continue
    ideal=math.dist([float(ads[a1][f'pdbx_model_Cartn_{c}_ideal']) for c in 'xyz'],[float(ads[a2][f'pdbx_model_Cartn_{c}_ideal']) for c in 'xyz'])
    actual=dist(lookup[a1],lookup[a2]); vals.append({'atoms':[a1,a2],'order':b['value_order'],'actual_A':actual,'ideal_A':ideal,'difference_A':actual-ideal})
   r['coordinate_bond_checks'].append({'structure':acc,'chain':lr['author_chain'],'bond_count':len(vals),'max_absolute_difference_A':max(abs(v['difference_A']) for v in vals),'bonds':vals})
 report['chemical_components'][ligand]=r
# Same-number sequence identity controls the mapping; inspect differences rather than silently matching different residues.
def ca(acc,ch):return {int(a['auth_seq_id']):a for a in structures[acc][1] if a['auth_asym_id']==ch and a['auth_atom_id']=='CA' and a['group_PDB']=='ATOM' and a['label_alt_id'] in ('.','A')}
def align(a,b):
 common=sorted(set(a)&set(b)); fit=[i for i in common if a[i]['auth_comp_id']==b[i]['auth_comp_id'] and i not in range(9,20)]
 x=np.array([a[i]['xyz'] for i in fit]);y=np.array([b[i]['xyz'] for i in fit]); cx=x.mean(0);cy=y.mean(0);u,s,vt=np.linalg.svd((x-cx).T@(y-cy)); rotation=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
 transformed=(x-cx)@rotation+cy; rmsd=float(np.sqrt(np.mean(np.sum((transformed-y)**2,axis=1))))
 assert np.allclose(rotation.T@rotation,np.eye(3),atol=1e-10) and np.linalg.det(rotation)>0.999
 return {'fit_CA_count':len(fit),'fit_CA_RMSD_A':rmsd,'sequence_mismatches':[{'residue':i,'first':a[i]['auth_comp_id'],'second':b[i]['auth_comp_id']} for i in common if a[i]['auth_comp_id']!=b[i]['auth_comp_id']], 'loop_CA_displacement_A':{i:float(np.linalg.norm((a[i]['xyz']-cx)@rotation+cy-b[i]['xyz'])) for i in common if i in range(9,20)}},rotation,cx,cy
for ac,ch,bc,dh in [('5NEV','A','5NEV','C'),('5LQF','A','5LQF','D'),('5NEV','A','5LQF','A'),('5NEV','C','5LQF','D')]:
 r,rot,cx,cy=align(ca(ac,ch),ca(bc,dh));r['first']=f'{ac}:{ch}';r['second']=f'{bc}:{dh}'
 for atom in ('CA','OH','CZ'):
  aa=[a for a in structures[ac][1] if a['auth_asym_id']==ch and a['auth_seq_id']=='15' and a['auth_atom_id']==atom]
  bb=[a for a in structures[bc][1] if a['auth_asym_id']==dh and a['auth_seq_id']=='15' and a['auth_atom_id']==atom]
  r['Tyr15_'+atom+'_displacement_A']=None if not aa or not bb else float(np.linalg.norm((aa[0]['xyz']-cx)@rot+cy-bb[0]['xyz']))
 report['comparisons'].append(r)
report['runtime']={'python':platform.python_version(),'biopython':Bio.__version__,'numpy':np.__version__,'elapsed_seconds':time.time()-start,'cpu_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,'maxrss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'threads':1,'device':'CPU','stop':'normal exit'}
(out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({'comparisons':report['comparisons'],'runtime':report['runtime']},indent=2))
