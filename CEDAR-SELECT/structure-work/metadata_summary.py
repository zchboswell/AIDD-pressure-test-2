import json
r=json.load(open('CEDAR-SELECT/structure-work/analysis.json'))
for a,s in r['structures'].items():
 print(a,'revision',s['revision'],'assembly',s['assemblies'])
 print('site missing residues',[x for x in s['unobserved_residues'] if x.get('auth_asym_id') in ('A','C' if a=='5NEV' else 'D') and int(x['auth_seq_id']) in list(range(9,20))+[31,33,64,80,81,82,83,84,85,86,89,131,132,134,135,145,146]])
 print('site missing atoms',[x for x in s['unobserved_atoms'] if x.get('auth_asym_id') in ('A','C' if a=='5NEV' else 'D') and int(x['auth_seq_id']) in list(range(9,20))+[31,33,64,80,81,82,83,84,85,86,89,131,132,134,135,145,146]])
 print('polymer spans',[(p['pdbx_strand_id'],len(p['pdbx_seq_one_letter_code_can'].replace('\n',''))) for p in s['polymers']])
for a,c in r['chemical_components'].items():print(a,[(v['chain'],v['max_absolute_difference_A']) for v in c['coordinate_bond_checks']])
