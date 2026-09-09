import json
p=json.load(open('CEDAR-SELECT/inputs/coverage/analogue_evidence.json'))
for c in p['compounds']:
 print(c['public_source_identifiers']['example'],c['smiles'])
 print([(o['target'],o['original_value']) for o in p['observations'] if o['compound_id']==c['compound_id']])
q=json.load(open('CEDAR-SELECT/inputs/public_assay_context.json'))
print(q.keys())
