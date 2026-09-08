"""Compare deterministic scientific outputs from the independent fresh-directory replay."""
import argparse,hashlib,json,pathlib
ROOT=pathlib.Path(__file__).resolve().parents[1]
a=argparse.ArgumentParser();a.add_argument('--replay',required=True,type=pathlib.Path);r=a.parse_args().replay
paths=['delivery/curation.csv','delivery/predictions.csv','delivery/selection.csv','delivery/experiments.csv','delivery/structure_checks.json','analysis/curation_details.csv','analysis/model_oof.csv','analysis/model_partitions.csv','analysis/model_fold_metrics.csv','analysis/model_summary.json','analysis/cycle_accounting.csv','analysis/cycle_summary.json','analysis/structure_evidence.json','analysis/structure_inventory_identity.json','analysis/legacy_split_overlap.json']
checks=[]
for p in paths:
 x=hashlib.sha256((ROOT/p).read_bytes()).hexdigest();y=hashlib.sha256((r/p).read_bytes()).hexdigest();checks.append({'path':p,'source_sha256':x,'replay_sha256':y,'identical':x==y})
result={'replay_root':str(r),'all_identical':all(c['identical'] for c in checks),'checks':checks,'scope':'Exact file equality of listed deterministic scientific tables/JSON; no claim of byte-identical environments or native capture.'}
(ROOT/'analysis/replay_comparison.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'all_identical':result['all_identical'],'files':len(checks)}));assert result['all_identical']
