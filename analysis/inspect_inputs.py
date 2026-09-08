import csv,json,collections
r=list(csv.DictReader(open('inputs/data/activity_archive.csv')))
for k in ['assay_id','endpoint','relation','unit','data_validity','record_origin']: print(k,dict(collections.Counter(x[k] for x in r)))
print('assay endpoints',dict(collections.Counter((x['assay_id'],x['endpoint'],x['relation']) for x in r)))
print('flagged rows',[(x['row_id'],x['data_validity']) for x in r if x['data_validity']])
print('primary values',[(x['row_id'],x['value'],x['unit']) for x in r if x['assay_id']=='CHEMBL5736732'][0:100])
refs=json.load(open('inputs/references/chembl_selected_records.json'));print('refs',len(refs));print('activity comments',[(x['activity_id'],x['activity_comment'],x['standard_value']) for x in refs[:8]])
