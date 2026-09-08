"""Reconcile all supplied rows without altering source bytes."""
import argparse,csv,hashlib,json,math,pathlib,collections
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from assay_math import normalize
ROOT=pathlib.Path(__file__).resolve().parents[1]
def readcsv(path):
 with open(path,newline='') as f:return list(csv.DictReader(f))
def writecsv(path,rows,fields=None):
 path.parent.mkdir(parents=True,exist_ok=True)
 with path.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)
def chemical(smiles):
 m=Chem.MolFromSmiles(smiles)
 return (m,Chem.MolToSmiles(m,isomericSmiles=True)) if m is not None else (None,None)
def digest(s):return hashlib.sha256(s.encode()).hexdigest()[:16]
def run(out):
 rows=readcsv(ROOT/'inputs/data/activity_archive.csv');refs={str(r['activity_id']):r for r in json.loads((ROOT/'inputs/references/chembl_selected_records.json').read_text())}
 manifest=json.loads((ROOT/'inputs/MANIFEST.json').read_text());assert all(hashlib.sha256((ROOT/'inputs'/p).read_bytes()).hexdigest()==sha for p,sha in manifest['files'].items())
 seen={};curated=[];sourcechecks=[];details=[];source_graphs=collections.defaultdict(set)
 for raw in refs.values():
  _,s=chemical(raw['canonical_smiles']);source_graphs[raw['molecule_chembl_id']].add(s)
 for r in rows:
  m,s=chemical(r['smiles']);n=normalize(r['value'],r['unit'],r['relation'],r['endpoint']);reason=[];status='context_only';sourceok=True
  if r['source_activity_id']:
   ref=refs.get(r['source_activity_id']);issues=[]
   if ref is None:issues.append('missing original source record')
   else:
    _,rs=chemical(ref['canonical_smiles']);rn=normalize(ref['standard_value'],ref['standard_units'],ref['standard_relation'],ref['standard_type'])
    for key,expected in [('chemical_id',ref['molecule_chembl_id']),('assay_id',ref['assay_chembl_id']),('endpoint',ref['standard_type']),('source_document',ref['document_chembl_id'])]:
     if r[key]!=str(expected):issues.append(key+' mismatch')
    if s!=rs:issues.append('chemical graph mismatch')
    if n['status']!='ok' or rn['status']!='ok' or not math.isclose(n.get('value_nm',0),rn.get('value_nm',1),rel_tol=1e-10) or n.get('relation_nm')!=rn.get('relation_nm'):issues.append('normalized source quantity mismatch')
   sourceok=not issues;sourcechecks.append({'row_id':r['row_id'],'source_activity_id':r['source_activity_id'],'source_matches':sourceok,'issues':'; '.join(issues)})
   if issues:reason+=issues
  if m is None:status='quarantined';reason.append('invalid SMILES; no supported chemical entity')
  elif n['status']!='ok':status='quarantined';reason.append(n['reason'])
  elif not sourceok:status='quarantined'
  elif r['data_validity']:status='quarantined';reason.append('unresolved source validity flag: '+r['data_validity'])
  elif r['source_activity_id'] and r['source_activity_id'] in seen:
   status='duplicate';reason.append('same source_activity_id as '+seen[r['source_activity_id']])
  elif r['assay_id']=='CHEMBL5736732' and r['endpoint']=='IC50' and n['relation_nm']=='=' and r['source_activity_id']:
   status='eligible';reason.append('source-reconciled exact historical assay-export IC50; not cellular or matched 1mM ATP endpoint')
  else:
   if r['assay_id']=='SIM-QC-CHECK':reason.append('simulated export example; no qualified biological assay; not training evidence')
   elif r['assay_id']!='CHEMBL5736732':reason.append('distinct assay context; no pooling with declared historical modeling endpoint')
   if n['relation_nm']!='=':reason.append('censored concentration; log threshold retained with reversed inequality, not exact training label')
   if r['endpoint']!='IC50':reason.append('distinct '+r['endpoint']+' endpoint')
  if status not in ('duplicate','quarantined') and r['source_activity_id']:seen[r['source_activity_id']]=r['row_id']
  if r['potential_duplicate']=='1':reason.append('source potential_duplicate flag retained; does not itself establish ingestion duplication or independence')
  valid=n['status']=='ok' and m is not None and sourceok and status!='quarantined'
  c={'row_id':r['row_id'],'status':status,'endpoint':r['endpoint'],'value_nm':n['value_nm'] if valid else '', 'relation_nm':n['relation_nm'] if valid else '', 'p_activity':f"{n['p_activity']:.9f}" if valid else '', 'p_relation':n['p_relation'] if valid else '', 'chemical_group':'CHEM-'+digest(s) if s else '', 'assay_group':r['assay_id']+'|'+r['endpoint'],'reason':'; '.join(reason)}
  curated.append(c)
  details.append({**c,'chemical_id':r['chemical_id'],'canonical_isomeric_smiles':s or '', 'murcko_scaffold':Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(m),isomericSmiles=True) if m else '', 'source_activity_id':r['source_activity_id'],'source_document':r['source_document'],'source_year':r['source_year'],'raw_value':r['value'],'raw_unit':r['unit'],'raw_relation':r['relation'],'raw_origin':r['record_origin'],'potential_duplicate':r['potential_duplicate']})
 inv=readcsv(ROOT/'inputs/data/inventory.csv');readouts={r['compound_id']:r for r in readcsv(ROOT/'inputs/data/program_readouts.csv')};vials={};canonical_by_id={};inventory=[];predictions=[]
 for r in inv:
  m,s=chemical(r['smiles']);status='valid';reason=[]
  if not s:status='quarantined';reason.append('invalid SMILES')
  elif r['source_id'] in source_graphs and s not in source_graphs[r['source_id']]:status='quarantined';reason.append('graph conflicts with source_id chemical identity')
  elif r['compound_id'] in canonical_by_id and s!=canonical_by_id[r['compound_id']]:status='quarantined';reason.append('graph conflicts with earlier compound_id representation')
  elif r['vial_id'] in vials:
   representative=vials[r['vial_id']]
   if s!=representative['canonical_isomeric_smiles']:status='quarantined';reason.append('vial identity conflict')
   else:status='duplicate';reason.append('same physical vial as '+representative['request_id']+'; no extra stock')
  else:reason.append('valid graph and no unresolved source-ID conflict')
  entry={**r,'identity_status':status,'canonical_isomeric_smiles':s or '', 'chemical_group':'CHEM-'+digest(s) if s else '', 'identity_reason':'; '.join(reason)}
  inventory.append(entry)
  if status=='valid':vials[r['vial_id']]=entry;canonical_by_id[r['compound_id']]=s
  p={'request_id':r['request_id'],'status':status if status in ('quarantined','duplicate') else 'abstained','endpoint':'cellular_CDK2_target_engagement_IC50','assay_id':'SIM-TE-01','p_activity':'','lower':'','upper':'','reason':'; '.join(reason)}
  ro=readouts.get(r['compound_id'],{})
  if status=='valid' and ro.get('cell_target_engagement_ic50_nM'):
   p['status']='observed';p['p_activity']=f"{9-math.log10(float(ro['cell_target_engagement_ic50_nM'])):.9f}";p['reason']='Supplied simulated SIM-TE-01 measurement; no replicate uncertainty supplied. Not a published observation or model prediction.'
  elif status=='valid':p['reason']='SIM-TE-01 unmeasured; no validated cellular model or transferable biochemical calibration. Measure ENGAGE if selected.'
  predictions.append(p)
 writecsv(out/'delivery/curation.csv',curated);writecsv(out/'analysis/curation_details.csv',details);writecsv(out/'analysis/source_reconciliation.csv',sourcechecks);writecsv(out/'analysis/inventory_identity.csv',inventory);writecsv(out/'delivery/predictions.csv',predictions)
 summary={'activity_rows':len(rows),'statuses':dict(collections.Counter(r['status'] for r in curated)),'eligible_unique_chemical_groups':len({r['chemical_group'] for r in curated if r['status']=='eligible'}),'source_comparisons':len(sourcechecks),'source_mismatches':[r for r in sourcechecks if not r['source_matches']],'inventory_statuses':dict(collections.Counter(r['identity_status'] for r in inventory)),'prediction_statuses':dict(collections.Counter(r['status'] for r in predictions)),'manifest_verified':True,'normalization':'preserve full canonical isomeric graph; no salt/tautomer/stereo/isotope stripping','model_endpoint':'historical CHEMBL5736732 exact IC50 export; cyclin/ATP per-row unresolved; not current cellular assay'}
 (out/'analysis/curation_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--out',type=pathlib.Path,default=ROOT);args=ap.parse_args();run(args.out.resolve())
