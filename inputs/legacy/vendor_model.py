"""Simulated notebook shortcut. CSV output and printed score are unvalidated."""
import csv, json, pathlib
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score, mean_absolute_error
from assay_math import normalize

def main(root):
    root=pathlib.Path(root)
    activity={r['row_id']:r for r in csv.DictReader((root/'data/activity_archive.csv').open())}
    x=[];y=[]
    for f in csv.DictReader((root/'data/vendor_features.csv').open()):
        r=activity[f['row_id']]
        if not f['vendor_activity_index']:continue
        n=normalize(r['value'],r['unit'],r['relation'],r['endpoint'])
        if n['status']!='ok':continue
        x.append([float(f[k]) for k in ['molecular_weight','logp','heavy_atoms','vendor_activity_index']]);y.append(n['p_activity'])
    a,b,c,d=train_test_split(np.asarray(x),np.asarray(y),test_size=.25,random_state=7)
    m=RandomForestRegressor(n_estimators=80,random_state=7,n_jobs=1).fit(a,c)
    p=m.predict(b);print(json.dumps({'r2':r2_score(d,p),'mae':mean_absolute_error(d,p),'rows':len(y),'split':'random_rows'}))

if __name__=='__main__':
    import sys
    main(sys.argv[1])
