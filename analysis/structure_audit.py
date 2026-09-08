"""Immutable CEDAR export/frame audit. Run --out NEW_DIRECTORY for replay."""
import argparse, csv, hashlib, json, math, pathlib, sys
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolAlign
P=pathlib.Path(__file__).resolve().parents[1]
def load(p):
    m=Chem.MolFromMolFile(str(p),removeHs=False)
    if m is None: raise ValueError(f'parse failed: {p}')
    assert all(a.GetAtomicNum()>1 for a in m.GetAtoms())
    assert np.isfinite(m.GetConformer().GetPositions()).all()
    return m

def maps(prb,ref):
    if prb.GetNumAtoms()!=ref.GetNumAtoms() or prb.GetNumBonds()!=ref.GetNumBonds(): raise ValueError('incomplete graph')
    def sig(a): return (a.GetAtomicNum(),a.GetIsotope(),a.GetFormalCharge(),a.GetTotalNumHs(),a.GetIsAromatic())
    matches=ref.GetSubstructMatches(prb,uniquify=False,useChirality=True,maxMatches=100000)
    if len(matches)>=100000: raise ValueError('enumeration limit')
    valid=[list(enumerate(x)) for x in matches if all(sig(prb.GetAtomWithIdx(i))==sig(ref.GetAtomWithIdx(j)) for i,j in enumerate(x))]
    if not valid: raise ValueError('no state-preserving bijection')
    return valid

def rms(prb,ref):
    mm=maps(prb,ref); x=prb.GetConformer().GetPositions(); y=ref.GetConformer().GetPositions()
    vals=[float(np.sqrt(np.mean([np.sum((x[i]-y[j])**2) for i,j in m]))) for m in mm]
    k=int(np.argmin(vals)); api=rdMolAlign.CalcRMS(prb,ref,map=mm,symmetrizeConjugatedTerminalGroups=False)
    assert abs(api-vals[k])<1e-9
    return dict(rmsd_A=api,map=mm[k],map_count=len(mm),max_mapped_displacement_A=float(max(np.linalg.norm(x[i]-y[j]) for i,j in mm[k])))

def raw(p):
    lines=p.read_text().splitlines(); n=int(lines[3][:3]); return [tuple(float(l[k:k+10]) for k in (0,10,20)) for l in lines[4:4+n]]
def scalar(x,y,m): return math.sqrt(sum(sum((x[i][k]-y[j][k])**2 for k in range(3)) for i,j in m)/len(m))
def transform(m,R,t):
    q=Chem.Mol(m)
    for i,p in enumerate(m.GetConformer().GetPositions()): q.GetConformer().SetAtomPosition(i,tuple(R@p+t))
    return q

def main():
    a=argparse.ArgumentParser(); a.add_argument('--out',default=str(P)); args=a.parse_args(); out=pathlib.Path(args.out); (out/'analysis').mkdir(parents=True,exist_ok=True); (out/'delivery').mkdir(exist_ok=True)
    d=P/'inputs/structures'; files=[d/x for x in ['1h1q_2a6_A.sdf','vendor_reference_export.sdf','vendor_pose_priority.sdf']]; ref,export,priority=[load(p) for p in files]
    results={k:rms(m,ref) for k,m in [('reference_export',export),('priority_export',priority)]}
    for key,p in zip(results,files[1:]):
        r=results[key]; r['independent_raw_scalar_rmsd_A']=scalar(raw(p),raw(files[0]),r['map']); assert abs(r['rmsd_A']-r['independent_raw_scalar_rmsd_A'])<1e-9
    R=np.array([[0,-1,0],[1,0,0],[0,0,1]]); t=np.array([12.,-7.,3.]); shift=transform(ref,np.eye(3),np.array([5.,0,0]))
    controls={'unchanged':rms(ref,ref)['rmsd_A'],'renumbered':rms(Chem.RenumberAtoms(ref,list(range(23,-1,-1))),ref)['rmsd_A'],'ligand_only_5A_shift':rms(shift,ref)['rmsd_A'],'common_rigid_transform_priority':rms(transform(priority,R,t),transform(ref,R,t))['rmsd_A']}
    assert controls['unchanged']<1e-9 and controls['renumbered']<1e-9 and abs(controls['ligand_only_5A_shift']-5)<1e-9
    assert abs(controls['common_rigid_transform_priority']-results['priority_export']['rmsd_A'])<1e-9
    for key,m in [('changed_formal_charge',Chem.Mol(ref)),('missing_atom',Chem.RWMol(ref))]:
        if key=='missing_atom': m.RemoveAtom(0)
        else: m.GetAtomWithIdx(13).SetFormalCharge(1)
        try: maps(m,ref)
        except ValueError as e: controls[key]={'status':'expected rejection','reason':str(e)}
        else: raise AssertionError(key)
    # A distinct chiral synthetic control checks stereochemical matching, not CEDAR evidence.
    s1=Chem.MolFromSmiles('F[C@](Cl)(Br)I'); s2=Chem.MolFromSmiles('F[C@@](Cl)(Br)I')
    try: maps(s1,s2)
    except ValueError: controls['changed_stereo']={'status':'expected rejection'}
    else: raise AssertionError('stereo')
    records=[]
    for l in (d/'1h1q.pdb').read_text().splitlines():
        if l[:6].strip() in ('ATOM','HETATM'):
            records.append(dict(serial=int(l[6:11]),name=l[12:16].strip(),alt=l[16].strip(),res=l[17:20],chain=l[21],resi=l[22:27].strip(),xyz=[float(l[k:k+8]) for k in (30,38,46)],occupancy=float(l[54:60]),bfactor=float(l[60:66]),element=l[76:78].strip()))
    lig=[r for r in records if r['res']=='2A6' and r['chain']=='A']; assert len(lig)==24
    refxyz=raw(files[0]); pdbmap=[]
    for i,x in enumerate(refxyz):
        hits=[r for r in lig if sum((x[k]-r['xyz'][k])**2 for k in range(3))<1e-12 and r['element']==ref.GetAtomWithIdx(i).GetSymbol()]; assert len(hits)==1
        pdbmap.append({'sdf_index':i,**hits[0]})
    assert len({r['serial'] for r in pdbmap})==24
    rec=[r for r in records if r['chain'] in ('A','B') and r['res']!='2A6']
    contacts=[]; minimums={}; closest={}
    for label,m in [('deposited',ref),('priority',priority)]:
        xyz=m.GetConformer().GetPositions(); minimums[label]=min(float(np.linalg.norm(x-np.array(r['xyz']))) for x in xyz for r in rec if r['res']!='HOH')
        closest[label]=min((dict(ligand_index=i,partner=f"{r['chain']}:{r['res']}{r['resi']}:{r['name']}",distance_A=float(np.linalg.norm(x-np.array(r['xyz'])))) for i,x in enumerate(xyz) for r in rec if r['res']!='HOH'),key=lambda z:z['distance_A'])
        for i,x in enumerate(xyz):
            if m.GetAtomWithIdx(i).GetSymbol() not in ('N','O'):continue
            for r in rec:
                dist=float(np.linalg.norm(x-np.array(r['xyz'])))
                if r['element'] in ('N','O') and dist<=3.5:
                    contacts.append({'pose':label,'ligand_index':i,'ligand_atom':pdbmap[i]['name'] if label=='deposited' else 'export_'+str(i),'partner':f"{r['chain']}:{r['res']}{r['resi']}:{r['name']}",'distance_A':dist,'interpretation':'polar proximity only; no H-bond energetic claim'})
    with (out/'analysis/structure_contacts.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=['pose','ligand_index','ligand_atom','partner','distance_A','interpretation']);w.writeheader();w.writerows(contacts)
    # Report fitted shape only on a copy, after immutable primary checks.
    fit=rdMolAlign.GetBestRMS(Chem.Mol(priority),ref,map=maps(priority,ref),symmetrizeConjugatedTerminalGroups=False)
    mp=results['priority_export']['map']; deltas=np.array([priority.GetConformer().GetPositions()[i]-ref.GetConformer().GetPositions()[j] for i,j in mp])
    bondchanges=[]
    for b in priority.GetBonds():
        i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx(); mapp=dict(mp);xp=priority.GetConformer().GetPositions();xr=ref.GetConformer().GetPositions(); bondchanges.append(abs(float(np.linalg.norm(xp[i]-xp[j])-np.linalg.norm(xr[mapp[i]]-xr[mapp[j]]))))
    entry=json.loads((d/'1h1q_entry.json').read_text())
    evidence={'runtime':{'python':sys.version,'rdkit_version':rdkit.__version__,'rdkit_path':rdkit.__file__},'sources':{str(p.relative_to(P)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files+[d/'1h1q.pdb',d/'1h1q.cif',d/'1h1q_entry.json']},'export_checks':results,'controls':controls,'reference_pdb_atom_mapping':pdbmap,'priority_translation_mean_A':deltas.mean(axis=0).tolist(),'priority_translation_residual_max_A':float(np.abs(deltas-deltas.mean(axis=0)).max()),'priority_mapped_bond_length_difference_max_A':max(bondchanges),'priority_ligand_fitted_shape_rmsd_A_SECONDARY':fit,'minimum_ligand_protein_heavy_atom_distance_A':minimums,'closest_protein_pairs':closest,'ligand_occupancy_range':[min(r['occupancy'] for r in lig),max(r['occupancy'] for r in lig)],'ligand_bfactor_range_A2':[min(r['bfactor'] for r in lig),max(r['bfactor'] for r in lig)],'ligand_altlocs':sorted(set(r['alt'] for r in lig)),'entry_title':entry['struct']['title'],'revision':entry['rcsb_accession_info'],'reference_pdb_match':'24 unique element-matched identical coordinates in deposited chain A 2A6','assessment_scope':'No docking, coordinate repair, full plausibility grading or electron-density validation; original frame audited.','contacts_path':'analysis/structure_contacts.csv'}
    (out/'analysis/structure_evidence.json').write_text(json.dumps(evidence,indent=2)+'\n')
    (out/'analysis/structure_api.txt').write_text(f'RDKit {rdkit.__version__}\n'+rdMolAlign.CalcRMS.__doc__)
    delivery={'source_accession':'1H1Q','ligand_identity_chain':{'component':'2A6','name':'NU6094, 2-anilino-6-cyclohexylmethoxypurine','author_chain':'A','label_asym_id':'E','author_residue':'1301','biological_assembly':'1; chains A/B; identity operator','molecular_graph_smiles':Chem.MolToSmiles(ref)},'atom_mapping_method':'Complete state-preserving heavy-atom graph isomorphisms, symmetry-aware minimum; no alignment; 24 heavy atoms; RDKit CalcRMS with explicit validated maps, independent raw-coordinate scalar sums. PDB correspondence independently element/coordinate checked.','reference_export_fixed_frame_rmsd_A':results['reference_export']['rmsd_A'],'priority_export_fixed_frame_rmsd_A':results['priority_export']['rmsd_A'],'priority_usable_for_contacts':False,'priority_interpretation':'Whole ligand translated relative to unchanged receptor; supplied export cannot support deposited binding-site contacts. Fitting would erase the defect.','evidence_paths':['analysis/structure_evidence.json','analysis/structure_contacts.csv','analysis/structure_audit.py','analysis/structure_plan.md','analysis/structure_notes.md']}
    delivery['ligand_identity_chain']['author_residue']=lig[0]['resi']
    (out/'delivery/structure_checks.json').write_text(json.dumps(delivery,indent=2)+'\n')
    print(json.dumps({'results':results,'controls':controls,'minimums':minimums,'fit_secondary':fit,'translation':deltas.mean(axis=0).tolist(),'contacts':contacts},indent=2))
if __name__=='__main__':main()
