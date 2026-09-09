"""Bounded integrity/accounting checks; no new scientific calculations."""
from pathlib import Path
import json,hashlib,re
R=Path(__file__).resolve().parents[2];C=R/'CEDAR-SELECT';out={}
m=json.loads((R/'admin/R005_INPUT_MANIFEST.json').read_text())['files']
for rel,h in m.items():assert hashlib.sha256((C/'inputs'/rel).read_bytes()).hexdigest()==h,rel
out['immutable_input_hashes_checked']=len(m)
b=json.loads((C/'delivery/experiments.json').read_text());assert sum(x['cap_usd'] for x in b['budget'])==120000
w=b['budget'][0];assert sum(w[k] for k in ['paid_usd','reserved_usd','new_O2_usd','new_O1_usd','remaining_uncommitted_usd'])==10000
assert sum(b['nominations'][0]['broader_material_plan_mg'].values())==2
assert len(b['nominations'])==4 and sum(x['incremental_usd'] for x in b['coverage_schedule'])==3400
out['budget_material_nomination_checks']='PASS'
text=(C/'delivery/decision.md').read_text();out['decision_words']=len(text.split());assert out['decision_words']<=2000
broken=[]
for p in (C/'delivery').glob('*.md'):
 for dest in re.findall(r'\]\(([^)]+)\)',p.read_text()):
  if dest.startswith(('http:','https:','#')):continue
  q=p.parent/dest.split('#')[0]
  if not q.exists():broken.append({'file':str(p),'link':dest})
assert not broken,broken
out['delivery_links']='PASS'
r=json.loads((C/'operations/replay/replay_report.json').read_text());assert r['status']=='PASS' and r['all_subprocesses_exited']
s=json.loads((C/'model-work/results/summary.json').read_text());assert s['compounds']==60 and s['source_rows']==120 and not s['joint_model_led_selection']
d=json.loads((C/'model-work/results/design_diagnostics.json').read_text());assert all(x['same_stereochemical_core'] for x in d);assert all('abstain' in x['status'] for x in d[:2])
out['replay_and_reported_status']='PASS';out['status']='PASS'
(C/'delivery/package_validation.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
