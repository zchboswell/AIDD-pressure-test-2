"""Replay in an isolated temporary copy; never modify the frozen reference outputs."""
import os,sys,json,shutil,tempfile,subprocess,time,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; C=ROOT/'CEDAR-SELECT'; OUT=Path(sys.argv[1]).resolve() if len(sys.argv)>1 else Path(tempfile.mkdtemp(prefix='cedar2-replay-report-'))
OUT.mkdir(parents=True,exist_ok=True);runtime=os.environ.get('CEDAR_PYTHON','/home/zbos/.local/share/cedar-runtime/20260908/cedar-python');start=time.time()
work=Path(tempfile.mkdtemp(prefix='cedar2-isolated-'))
shutil.copytree(C,work/'CEDAR-SELECT',ignore=shutil.ignore_patterns('__pycache__','*.pyc','mpl-cache','MANIFEST.sha256'))
(work/'inputs/data').mkdir(parents=True);shutil.copy2(ROOT/'inputs/data/inventory.csv',work/'inputs/data/inventory.csv')
# Remove generated accepted references before execution; supplied source text derivatives remain inputs to structural review.
shutil.rmtree(work/'CEDAR-SELECT/model-work/results')
for name in ['coverage-work/analogue_reconciled.json','coverage-work/budget_material_schedule.json','coverage-work/local_pairs.json','coverage-work/retrieval_audit.json','coverage-work/model_validation.json','design-work/designs.json','structure-work/analysis.json','structure-work/sensitivity.json','structure-work/molecular_view.png']:
 (work/'CEDAR-SELECT'/name).unlink(missing_ok=True)
env=os.environ.copy();env.update({'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','MPLCONFIGDIR':str(work/'mpl-cache')});runs=[]
scripts=['model-work/run.py','coverage-work/analyze.py','design-work/designs.py','model-work/score_designs.py','coverage-work/model_check.py','structure-work/analyze.py','structure-work/sensitivity.py','structure-work/render.py','model-work/figures.py']
for script in scripts:
 t=time.time();p=subprocess.run([runtime,'CEDAR-SELECT/'+script],cwd=work,env=env,text=True,capture_output=True,timeout=600);(OUT/(script.replace('/','_')+'.log')).write_text(p.stdout+p.stderr);runs.append({'script':script,'exit_code':p.returncode,'wall_seconds':time.time()-t});assert p.returncode==0,(script,p.stderr[-2000:])
# Replay comparison excludes runtime/environment timing values only, not scientific arrays.
def strip(x):
 if isinstance(x,dict):return {k:strip(v) for k,v in x.items() if k not in ['runtime']}
 if isinstance(x,list):return [strip(v) for v in x]
 return x
checks=[]
paths=['model-work/results/summary.json','model-work/results/design_diagnostics.json','coverage-work/analogue_reconciled.json','coverage-work/budget_material_schedule.json','coverage-work/local_pairs.json','coverage-work/retrieval_audit.json','design-work/designs.json','structure-work/analysis.json','structure-work/sensitivity.json']
for rel in paths:
 a=strip(json.loads((C/rel).read_text()));b=strip(json.loads((work/'CEDAR-SELECT'/rel).read_text()));assert a==b,rel;checks.append({'path':rel,'scientific_values_match':True})
for p in sorted((C/'model-work/results').glob('*.csv')):
 assert p.read_bytes()==(work/'CEDAR-SELECT/model-work/results'/p.name).read_bytes(),p.name;checks.append({'path':str(p.relative_to(C)),'bytes_match':True})
for p in (C/'inputs').rglob('*'):
 if p.is_file():assert p.read_bytes()==(work/'CEDAR-SELECT'/p.relative_to(C)).read_bytes()
report={'status':'PASS','scripts':runs,'comparisons':checks,'input_bytes_unchanged':True,'elapsed_seconds':time.time()-start,'threads_per_job':1,'jobs_sequential':True,'all_subprocesses_exited':True,'temporary_copy':str(work),'temporary_copy_retained':False,'runtime':runtime,'replay_script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
shutil.rmtree(work);(OUT/'replay_report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
