"""Independent contract/physical-quantity/accounting checks on written artifacts."""
import argparse,collections,csv,decimal,hashlib,json,math,pathlib,re
ROOT=pathlib.Path(__file__).resolve().parents[1]
D=decimal.Decimal;decimal.getcontext().prec=40
def rows(p):
 with p.open(newline='') as f:return list(csv.DictReader(f))
def run(out):
 checks=[]
 def check(name,value):
  checks.append({'check':name,'pass':bool(value)})
  if not value:raise AssertionError(name)
 manifest=json.loads((ROOT/'inputs/MANIFEST.json').read_text());check('immutable input hashes',all(hashlib.sha256((ROOT/'inputs'/p).read_bytes()).hexdigest()==h for p,h in manifest['files'].items()))
 stock=rows(ROOT/'inputs/data/inventory.csv');activity=rows(ROOT/'inputs/data/activity_archive.csv');cur=rows(out/'delivery/curation.csv');pred=rows(out/'delivery/predictions.csv');sel=rows(out/'delivery/selection.csv');exp=rows(out/'delivery/experiments.csv');menu={r['experiment_id']:r for r in rows(ROOT/'inputs/data/experiment_menu.csv')};readouts={r['compound_id']:r for r in rows(ROOT/'inputs/data/program_readouts.csv')}
 def ids(table,key):return collections.Counter(r[key] for r in table)
 check('activity rows exactly once',ids(activity,'row_id')==ids(cur,'row_id') and all(v==1 for v in ids(cur,'row_id').values()))
 check('stock request accounting exactly once',ids(stock,'request_id')==ids(pred,'request_id')==ids(sel,'request_id'))
 check('curation tokens',all(r['status'] in {'eligible','context_only','duplicate','quarantined'} for r in cur));check('prediction tokens',all(r['status'] in {'predicted','observed','abstained','duplicate','quarantined'} for r in pred));check('selection tokens',all(r['status'] in {'advance','diagnostic','reserve','duplicate','quarantined'} for r in sel))
 unitmol={'nM':D('1e-9'),'uM':D('1e-6'),'µM':D('1e-6'),'μM':D('1e-6'),'mM':D('1e-3'),'M':D(1)};inverse={'=':'=','<':'>','<=':'>=','>':'<','>=':'<=','≤':'>=','≥':'<='};raw={r['row_id']:r for r in activity};quantities=0
 for r in cur:
  for key in ['value_nm','p_activity']:
   if r[key]:check(r['row_id']+' finite '+key,math.isfinite(float(r[key])))
  if r['p_activity']:
   source=raw[r['row_id']];mol=D(source['value'])*unitmol[source['unit']];p=-mol.log10();check(r['row_id']+' independent Decimal molarity/log oracle',abs(D(r['p_activity'])-p)<D('0.000000001') and abs(D(r['value_nm'])-mol*D('1e9'))<D('0.00001'));check(r['row_id']+' reversed relation',r['p_relation']==inverse[source['relation']]);check(r['row_id']+' endpoint intact',r['endpoint']==source['endpoint']);quantities+=1
  if r['status']=='eligible':check(r['row_id']+' eligibility scope',r['relation_nm']=='=' and raw[r['row_id']]['assay_id']=='CHEMBL5736732' and r['endpoint']=='IC50')
 st={r['request_id']:r for r in stock}
 for r in pred:
  check(r['request_id']+' prediction endpoint attribution',r['endpoint']=='cellular_CDK2_target_engagement_IC50' and r['assay_id']=='SIM-TE-01')
  if r['status'] in ('predicted','observed'):check(r['request_id']+' estimate or bound exists',any(r[k] for k in ['p_activity','lower','upper']))
  if r['status']=='observed':
   source=readouts[st[r['request_id']]['compound_id']]['cell_target_engagement_ic50_nM'];check(r['request_id']+' observed scalar source',abs(D(r['p_activity'])+ (D(source)*D('1e-9')).log10())<D('0.000000001'))
  for k in ['p_activity','lower','upper']:
   if r[k]:check(r['request_id']+' finite prediction '+k,math.isfinite(float(r[k])))
  if r['lower'] and r['upper']:check(r['request_id']+' ordered bounds',float(r['lower'])<=float(r['upper']))
 selected=[r for r in sel if r['status'] in ('advance','diagnostic')];check('roles at most6plus2',sum(r['status']=='advance' for r in selected)<=6 and sum(r['status']=='diagnostic' for r in selected)<=2);check('unique compounds and physical vials <=8',len({r['compound_id'] for r in selected})==len({st[r['request_id']]['vial_id'] for r in selected})==len(selected)<=8)
 check('every experiment refers to selected compound',all(e['compound_id'] in {r['compound_id'] for r in selected} for e in exp));check('one row per compound experiment',len({(r['compound_id'],r['experiment_id']) for r in exp})==len(exp));cost=D(0);account=[]
 for s in selected:
  r=st[s['request_id']];es=[e for e in exp if e['compound_id']==s['compound_id']];mass=D(0);local=D(r['handling_usd']);identity=0;duration=0
  for e in es:
   n=D(e['replicates']);check(s['compound_id']+' valid '+e['experiment_id'],e['experiment_id'] in menu and n>0 and n==n.to_integral_value() and bool(e['purpose']));m=menu[e['experiment_id']];mass+=D(m['mass_mg'])*n;local+=D(m['cost_usd'])*n
   if e['experiment_id']=='IDENTITY':identity+=int(m['turnaround_days'])*int(n)
   else:duration=max(duration,int(m['turnaround_days']))
  day=int(r['delivery_days'])+identity+duration;check(s['compound_id']+' mass/deadline/purity',mass<=D(r['mass_mg']) and day<=8 and (D(r['purity_pct'])>=90 or identity>=1));cost+=local;account.append({'compound_id':s['compound_id'],'cost_usd':str(local),'mass_mg':str(mass),'completion_day':day})
 check('budget <=18000',cost<=18000);saved=json.loads((out/'analysis/cycle_summary.json').read_text());check('saved cost agrees independent Decimal sum',cost==D(str(saved['total_usd'])))
 s=json.loads((out/'delivery/structure_checks.json').read_text());check('structure fields finite and original frame',math.isclose(float(s['reference_export_fixed_frame_rmsd_A']),0,abs_tol=1e-9) and math.isclose(float(s['priority_export_fixed_frame_rmsd_A']),math.sqrt(77),abs_tol=1e-9) and s['priority_usable_for_contacts'] is False);check('structure evidence exists',all((out/p).is_file() for p in s['evidence_paths']))
 required=['delivery/report.md','delivery/selection.csv','delivery/experiments.csv','delivery/curation.csv','delivery/predictions.csv','delivery/structure_checks.json','delivery/REPRODUCE.md','delivery/NEXT_BATCH.md']
 for p in required:check('required artifact '+p,(out/p).is_file() and (out/p).stat().st_size>0 and (out/p).stat().st_size<10_000_000)
 for name in ['decision_comparison.png','model_diagnostic.png','structure_comparison.png']:check('figure '+name,(out/'delivery/figures'/name).stat().st_size>10000)
 # Reporting links resolved against either output or immutable project input/code tree.
 for p in ['delivery/report.md','delivery/NEXT_BATCH.md']:
  for link in re.findall(r'\]\(([^)]+)\)',(out/p).read_text()):
   if '://' not in link:check('report link '+link,((out/p).parent/link).exists() or ((ROOT/p).parent/link).exists())
 result={'all_pass':True,'check_count':len(checks),'independent_normalized_rows':quantities,'total_cost_usd':str(cost),'compound_accounting':account,'scope':'Schema, independent Decimal normalization/accounting, frozen source hashes, numeric frame reference and artifact existence; not prospective biological validation.','checks':checks}
 (out/'analysis/validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:result[k] for k in ['all_pass','check_count','independent_normalized_rows','total_cost_usd']}))
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--out',type=pathlib.Path,default=ROOT);run(a.parse_args().out.resolve())
