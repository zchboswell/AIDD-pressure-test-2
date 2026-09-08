import json,pathlib,csv
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
root=pathlib.Path(__file__).resolve().parents[1]
ref=Chem.SDMolSupplier(str(root/'inputs/structures/1h1q_2a6_A.sdf'),removeHs=True)[0]
inv={r['compound_id']:r for r in csv.DictReader(open(root/'inputs/data/inventory.csv'))}
a=Chem.MolFromSmiles(inv['CED-23']['smiles']);taut=rdMolStandardize.TautomerEnumerator()
record={'reference_isomeric_smiles':Chem.MolToSmiles(ref),'CED23_isomeric_smiles':Chem.MolToSmiles(a),'exact_graph_state_matches':Chem.MolToSmiles(ref)==Chem.MolToSmiles(a),'canonical_tautomer_matches':Chem.MolToSmiles(taut.Canonicalize(ref))==Chem.MolToSmiles(taut.Canonicalize(a)),'interpretation':'CED23 shares a tautomer-normalized chemical parent with deposited ligand but is not the identical supplied atom-specific tautomer. Tautomer equivalence is a grouping diagnostic only, not proof of measured state or same pose. No inventory state or coordinate was altered.'}
(root/'analysis/structure_inventory_identity.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))
