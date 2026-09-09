"""Illustrative superposition of immutable deposited complexes, not a predicted pose."""
import os,time,resource,json
os.environ['MPLCONFIGDIR']='CEDAR-SELECT/structure-work/mpl-cache'
from pathlib import Path
import numpy as np
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SVDSuperimposer import SVDSuperimposer
import pymol
pymol.finish_launching(['pymol','-cq'])
from pymol import cmd
root=Path('CEDAR-SELECT/structure-work');start=time.time();r=json.loads((root/'analysis.json').read_text());c=r['comparisons'][2]
def ca(acc):
 d=MMCIF2Dict('CEDAR-SELECT/inputs/structural_sources/'+acc+'.cif'); keys=['auth_asym_id','auth_seq_id','auth_atom_id','auth_comp_id','Cartn_x','Cartn_y','Cartn_z']; vals=zip(*(d['_atom_site.'+k] for k in keys));return {int(i):(n,np.array([float(x),float(y),float(z)])) for ch,i,at,n,x,y,z in vals if ch=='A' and at=='CA'}
a=ca('5NEV');b=ca('5LQF');fit=[(int(i),j) for i,j in c['residue_mapping'].items() if 1<=int(i)<=90 and int(i) not in range(9,20) and a[int(i)][0]==b[j][0]];sv=SVDSuperimposer();sv.set(np.array([b[j][1] for i,j in fit]),np.array([a[i][1] for i,j in fit]));sv.run();rot,tran=sv.get_rotran()
cmd.load('CEDAR-SELECT/inputs/structural_sources/5NEV.cif','cdk2');cmd.load('CEDAR-SELECT/inputs/structural_sources/5LQF.cif','cdk1')
# CDK2 original atoms transformed for display only into CDK1 frame.
def trans(x,y,z):return tuple(np.array([x,y,z])@rot+tran)
cmd.alter_state(1,'cdk2','(x,y,z)=trans(x,y,z)',space={'trans':trans});cmd.hide('everything');cmd.bg_color('white');cmd.set('ray_opaque_background',1);cmd.set('orthoscopic',1);cmd.set('antialias',2);cmd.set('max_threads',1)
cmd.show('cartoon','(cdk1 or cdk2) and chain A and resi 1-90');cmd.color('gray85','cdk1');cmd.color('gray70','cdk2');cmd.color('marine','cdk2 and chain A and resi 9-19');cmd.color('orange','cdk1 and chain A and resi 9-19')
cmd.show('sticks','chain A and resi 15 and (cdk1 or cdk2)');cmd.set('stick_radius',0.18);cmd.show('sticks','cdk2 and resn 72L and chain A');cmd.color('forest','cdk2 and resn 72L and chain A');cmd.color('red','cdk2 and resn 72L and elem O');cmd.color('blue','cdk2 and resn 72L and elem N');cmd.color('yellow','cdk2 and resn 72L and elem S')
cmd.orient('cdk2 and chain A and (resn 72L or resi 9-19)');cmd.zoom('(cdk2 and chain A and (resn 72L or resi 9-19)) or (cdk1 and chain A and resi 15)',4);cmd.turn('x',20);cmd.png(str(root/'molecular_view.png'),1000,750,dpi=150,ray=1)
(root/'render_runtime.json').write_text(json.dumps({'wall_seconds':time.time()-start,'cpu_seconds':resource.getrusage(resource.RUSAGE_SELF).ru_utime+resource.getrusage(resource.RUSAGE_SELF).ru_stime,'maxrss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,'threads':1,'stop':'normal exit','pymol_version':cmd.get_version()[0],'display_transform_rotation':rot.tolist(),'display_transform_translation':tran.tolist()},indent=2)+'\n')
cmd.quit()
