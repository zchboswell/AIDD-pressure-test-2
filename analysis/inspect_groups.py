import csv,collections,json
r=[x for x in csv.DictReader(open('analysis/curation_details.csv')) if x['status']=='eligible'];groups=collections.defaultdict(list)
for x in r:groups[x['chemical_group']].append(x)
print('repeated chemistry',json.dumps([[{'row_id':x['row_id'],'chemical_id':x['chemical_id'],'p_activity':x['p_activity'],'source_activity_id':x['source_activity_id']} for x in v] for v in groups.values() if len(v)>1],indent=2))
print('scaffolds',len(set(x['murcko_scaffold'] for x in r)));print('group_sizes',sorted(collections.Counter(x['murcko_scaffold'] for x in r).values(),reverse=True))
