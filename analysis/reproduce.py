"""Replay into a new self-contained output tree; never overwrite source results."""
import argparse,datetime,hashlib,json,os,pathlib,shutil,subprocess,sys,time
ROOT=pathlib.Path(__file__).resolve().parents[1]
def run(out,skip):
 if out.exists():raise ValueError('Replay output must be a new directory: '+str(out))
 if ROOT==out or ROOT in out.parents:raise ValueError('Use a directory outside source repository to avoid recursive copies.')
 started=datetime.datetime.now(datetime.timezone.utc).isoformat();out.mkdir(parents=True);shutil.copytree(ROOT/'inputs',out/'inputs');(out/'analysis').mkdir();(out/'delivery').mkdir();(out/'analysis/checks').mkdir();(out/'analysis/replay_logs').mkdir()
 for p in (ROOT/'analysis').iterdir():
  if p.is_file() and p.suffix in {'.py','.md','.patch'}:shutil.copy2(p,out/'analysis'/p.name)
 for p in (ROOT/'delivery').glob('*.md'):shutil.copy2(p,out/'delivery'/p.name)
 launcher=os.environ.get('CEDAR_PYTHON','/home/zbos/.local/share/cedar-runtime/20260908/cedar-python')
 environment=os.environ.copy();environment.update({'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','MKL_NUM_THREADS':'1','NUMEXPR_NUM_THREADS':'1','MPLCONFIGDIR':str(out/'.cache/matplotlib')})
 steps=[];receipt={'started_utc':started,'output_root':str(out),'cpu_affinity':sorted(os.sched_getaffinity(0)),'threads':1,'molecular_render_regenerated':not skip,'steps':steps,'complete':False}
 def save():
  (out/'analysis/replay_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 def execute(name,args,expected=0,python=True):
  command=([launcher] if python else [])+args;t=time.monotonic();result=subprocess.run(command,cwd=out,env=environment,capture_output=True,text=True);log=out/'analysis/replay_logs'/f'{name}.log';log.write_text(result.stdout+'\n--- STDERR ---\n'+result.stderr);step={'name':name,'exit_code':result.returncode,'expected_exit_code':expected,'seconds':round(time.monotonic()-t,3),'log':str(log.relative_to(out))};steps.append(step);save();print(json.dumps(step),flush=True)
  if result.returncode!=expected:raise RuntimeError(name+' failed; see '+str(log))
 try:
  for script in ['curate','model_audit','legacy_split_audit','select_cycle','structure_audit','structural_inventory_identity']:
   execute(script,[f'analysis/{script}.py','--out',str(out)])
  execute('helper_baseline',['analysis/test_assay_math.py','--module','inputs/legacy/assay_math.py','--output','analysis/checks/assay_math_baseline.json'],expected=1)
  execute('helper_candidate',['analysis/test_assay_math.py','--module','analysis/assay_math.py','--output','analysis/checks/assay_math_candidate.json'])
  execute('helper_independent',['analysis/independent_assay_oracle.py','--out',str(out)])
  execute('patch_reconstruction',['patch','--output=analysis/checks/assay_math_reconstructed.py','inputs/legacy/assay_math.py','analysis/assay_math.patch'],python=False)
  assert (out/'analysis/assay_math.py').read_bytes()==(out/'analysis/checks/assay_math_reconstructed.py').read_bytes()
  receipt['patch_reconstruction_sha256']=hashlib.sha256((out/'analysis/checks/assay_math_reconstructed.py').read_bytes()).hexdigest()
  execute('decision_figures',['analysis/plot_delivery.py','--out',str(out)])
  if skip:
   for p in (ROOT/'delivery/figures').glob('structure*'):shutil.copy2(p,out/'delivery/figures'/p.name)
  else:execute('molecular_figures',['analysis/structure_render.py','--out',str(out)])
  execute('validation',['analysis/validate_delivery.py','--out',str(out)])
  receipt['complete']=True;receipt['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();save()
 except Exception as exc:
  receipt['error']=repr(exc);receipt['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();save();raise
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--out',required=True,type=pathlib.Path);a.add_argument('--skip-molecular-render',action='store_true');args=a.parse_args();run(args.out.resolve(),args.skip_molecular_render)
