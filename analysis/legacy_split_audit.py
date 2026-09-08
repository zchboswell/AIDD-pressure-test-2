"""Retain actual original row-split overlap; never use this split for deployment."""
import argparse,collections,importlib.util,json,pathlib
import numpy as np
from sklearn.model_selection import train_test_split
from rdkit import Chem
from curate import readcsv
ROOT=pathlib.Path(__file__).resolve().parents[1]
a=argparse.ArgumentParser();a.add_argument('--out',type=pathlib.Path,default=ROOT);out=a.parse_args().out
spec=importlib.util.spec_from_file_location('original_assay_helper',ROOT/'inputs/legacy/assay_math.py');legacy=importlib.util.module_from_spec(spec);spec.loader.exec_module(legacy)
activity={r['row_id']:r for r in readcsv(ROOT/'inputs/data/activity_archive.csv')};rows=[]
for f in readcsv(ROOT/'inputs/data/vendor_features.csv'):
 r=activity[f['row_id']]
 if f['vendor_activity_index'] and legacy.normalize(r['value'],r['unit'],r['relation'],r['endpoint'])['status']=='ok':
  m=Chem.MolFromSmiles(r['smiles']);rows.append({**r,'graph':Chem.MolToSmiles(m,isomericSmiles=True) if m else 'INVALID:'+r['row_id']})
train,test=train_test_split(np.arange(len(rows)),test_size=.25,random_state=7)
result={'rows':len(rows),'train_rows':len(train),'test_rows':len(test),'overlapping_source_activity_ids':sorted(({rows[i]['source_activity_id'] for i in train}&{rows[i]['source_activity_id'] for i in test})-{''}),'overlapping_chemical_graphs':len({rows[i]['graph'] for i in train}&{rows[i]['graph'] for i in test}),'test_rows_with_training_chemical':[rows[i]['row_id'] for i in test if rows[i]['graph'] in {rows[j]['graph'] for j in train}],'mixed_assay_counts':dict(collections.Counter(r['assay_id'] for r in rows)),'interpretation':'Actual supplied vendor random-row split. Overlap is not prospective chemical-group validation; feature provenance and endpoint mismatch independently block deployment.'}
(out/'analysis/legacy_split_overlap.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
