"""Root-authored analytic references independent of candidate implementation."""
import hashlib,json,math,pathlib,argparse
a=argparse.ArgumentParser();a.add_argument("--out",type=pathlib.Path,default=pathlib.Path(__file__).resolve().parents[1]);out=a.parse_args().out
from assay_math import normalize
# Exact molarity powers and independent base-10 reference ln(2)/ln(10).
cases=[('M power',('1e-8','M','=','Kd'),10.,8.,'='),('mM power',('1e-5','mM','>','IC50'),10.,8.,'<'),('micro bound',('0.1','μM','≤','Ki'),100.,7.,'>='),('nano bound',('1000','nM','>=','Kd'),1000.,6.,'<='),('nonpower',('200','nM','<','IC50'),200.,7.-math.log(2)/math.log(10),'>')]
checks=[]
for name,args,nm,p,rel in cases:
 r=normalize(*args);checks.append({'name':name,'observed':r,'pass':r['status']=='ok' and math.isclose(r['value_nm'],nm,rel_tol=1e-12) and math.isclose(r['p_activity'],p,abs_tol=1e-12) and r['p_relation']==rel and r['endpoint']==args[3]})
for value,unit,relation,endpoint in [('inf','nM','=','IC50'),('NaN','nM','=','Ki'),('0','M','=','Kd'),('-1','uM','=','IC50'),('1','ng/mL','=','IC50'),('1','nM','!=','IC50'),('1','nM','=','EC50')]:
 r=normalize(value,unit,relation,endpoint);checks.append({'name':str((value,unit,relation,endpoint)),'observed':r,'pass':r['status']=='invalid'})
p=pathlib.Path(__file__).with_name('assay_math.py')
result={'candidate_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'oracle':'independent root-authored analytic molarity powers, separately expressed log identity and explicit contract rejection cases; non-blind review','checks':checks,'all_pass':all(c['pass'] for c in checks)}
path=out/'analysis/checks/assay_independent_oracle.json';path.parent.mkdir(parents=True,exist_ok=True);path.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({'all_pass':result['all_pass'],'cases':len(checks),'path':str(path)}))
assert result['all_pass']
