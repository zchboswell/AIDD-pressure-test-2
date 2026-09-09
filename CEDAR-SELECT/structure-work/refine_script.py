from pathlib import Path
p=Path('CEDAR-SELECT/structure-work/analyze.py');s=p.read_text();s=s.replace('protein comparison on common same-number identical-residue CA excluding residues 9-19;','protein comparison on sequence-aligned identical-residue CA excluding residues 9-19;')
a=s.index('def align(a,b):');b=s.index("for ac,ch,bc,dh",a)
s=s[:a]+'''def align(a,b):
 from Bio.Align import PairwiseAligner, substitution_matrices
 from Bio.SeqUtils import seq1
 ai=sorted(a);bi=sorted(b)
 sa=''.join(seq1(a[i]['auth_comp_id'],custom_map={'TPO':'T'}) for i in ai)
 sb=''.join(seq1(b[i]['auth_comp_id'],custom_map={'TPO':'T'}) for i in bi)
 aligner=PairwiseAligner();aligner.substitution_matrix=substitution_matrices.load('BLOSUM62');aligner.open_gap_score=-10;aligner.extend_gap_score=-0.5
 al=aligner.align(sa,sb)[0];mapping={ai[int(i)]:bi[int(j)] for i,j in zip(*al.indices) if i>=0 and j>=0}
 fit=[(i,j) for i,j in mapping.items() if a[i]['auth_comp_id']==b[j]['auth_comp_id'] and i not in range(9,20)]
 x=np.array([a[i]['xyz'] for i,j in fit]);y=np.array([b[j]['xyz'] for i,j in fit]);cx=x.mean(0);cy=y.mean(0);u,sv,vt=np.linalg.svd((x-cx).T@(y-cy));rotation=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
 transformed=(x-cx)@rotation+cy;rmsd=float(np.sqrt(np.mean(np.sum((transformed-y)**2,axis=1))))
 assert np.allclose(rotation.T@rotation,np.eye(3),atol=1e-10) and np.linalg.det(rotation)>0.999
 # Independent Bio.PDB implementation verifies superposition numerics and orientation.
 from Bio.SVDSuperimposer import SVDSuperimposer
 checker=SVDSuperimposer();checker.set(y,x);checker.run();assert abs(checker.get_rms()-rmsd)<1e-10
 return {'alignment':str(al),'residue_mapping':mapping,'fit_CA_count':len(fit),'fit_CA_RMSD_A':rmsd,'sequence_mismatches':[{'first_residue':i,'second_residue':j,'first':a[i]['auth_comp_id'],'second':b[j]['auth_comp_id']} for i,j in mapping.items() if a[i]['auth_comp_id']!=b[j]['auth_comp_id']], 'loop_CA_displacement_A':{i:float(np.linalg.norm((a[i]['xyz']-cx)@rotation+cy-b[j]['xyz'])) for i,j in mapping.items() if i in range(9,20)}},rotation,cx,cy
''' +s[b:]
s=s.replace("print(json.dumps({'comparisons':report['comparisons'],'runtime':report['runtime']},indent=2))","print(json.dumps({'summary':[{'pair':r['first']+' '+r['second'],'n':r['fit_CA_count'],'rmsd':r['fit_CA_RMSD_A'],'Tyr15_OH':r['Tyr15_OH_displacement_A']} for r in report['comparisons']], 'runtime':report['runtime']},indent=2))")
p.write_text(s)
