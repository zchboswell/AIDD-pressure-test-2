import json
r=json.load(open('CEDAR-SELECT/structure-work/analysis.json'))
for k,s in r['structures'].items():
 print(k,'entities',[(x['id'],x['pdbx_description'],x.get('pdbx_mutation')) for x in s['entities']]); print('seqdiff',s['sequence_differences'])
 for l in s['ligands']:
  print('lig',l['author_chain'],l['label_chain'],l['heavy_atom_count'],l['occupancy_range']); print('contacts',l['residue_contact_minima']);print('polar',l['polar_proximities']);print('water',l['nearby_waters'])
for l in r['chemical_components'].values(): print(l['name'],l['descriptors'])
for c in r['comparisons'][2:]: print(c['alignment'])
