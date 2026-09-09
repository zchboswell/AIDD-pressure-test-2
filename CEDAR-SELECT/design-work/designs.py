import json,csv,time,resource
from rdkit import Chem
from rdkit.Chem import Descriptors,Crippen,rdMolDescriptors
start=time.time()
p=json.load(open('CEDAR-SELECT/inputs/coverage/analogue_evidence.json'))
q=json.load(open('CEDAR-SELECT/inputs/public_assay_context.json'))
rows=list(csv.DictReader(open('inputs/data/inventory.csv')))
known={}; invalid=[]
for c in p['compounds']+q['train']+q['development']+rows:
 m=Chem.MolFromSmiles(c['smiles'])
 if m is None: invalid.append(c['compound_id']); continue
 known.setdefault(Chem.MolToSmiles(m),[]).append(c['compound_id'])
base='O=C1N([C@@H]2CCC[C@@H](O)C2)c2nc(Nc3c[nH]nc3%s)ncc2C12CC2'
def record(id,sub,status):
 s=base%sub;m=Chem.MolFromSmiles(s);assert m
 return dict(design_id=id,smiles=s,canonical_smiles=Chem.MolToSmiles(m),role=status,supplied_matches=sorted(set(known.get(Chem.MolToSmiles(m),[]))),MW=Descriptors.MolWt(m),cLogP=Crippen.MolLogP(m),TPSA=rdMolDescriptors.CalcTPSA(m),HBD=rdMolDescriptors.CalcNumHBD(m),HBA=rdMolDescriptors.CalcNumHBA(m),stereocenters=Chem.FindMolChiralCenters(m,includeUnassigned=True))
out=[record('SEL-D01','OCC(F)(F)CO','unmeasured exploratory synthesis'),record('SEL-D02','OCC(F)(F)C','unmeasured exploratory synthesis'),record('SEL-C39','OCCF','historical Example39 diagnostic comparator'),record('SEL-C41','OCC(F)F','historical Example41 diagnostic comparator')]
json.dump(out,open('CEDAR-SELECT/design-work/designs.json','w'),indent=2)
json.dump(dict(job='identity/descriptors only',threads=1,elapsed_seconds=time.time()-start,maxrss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,finished=True,invalid_identity_rows=invalid,failed_attempts=['python unavailable in PATH','cedar-python does not support stdin hyphen; script execution used','helper inspect.py shadowed standard library; renamed','invalid inventory SMILES caused first identity pass failure; excluded explicitly']),open('CEDAR-SELECT/design-work/resource_log.json','w'),indent=2)
print(json.dumps(out,indent=2))
