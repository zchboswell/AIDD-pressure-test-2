import os
os.environ['MPLCONFIGDIR']='/tmp/cedar2-mpl'
import json,csv
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import Draw
R=Path(__file__).resolve().parents[1];out=R/'delivery/figures';out.mkdir(parents=True,exist_ok=True)
a=json.loads((R/'inputs/coverage/analogue_evidence.json').read_text());points=[]
for c in a['compounds']:
 obs={r['target']:1000*float(r['original_value']) for r in a['observations'] if r['compound_id']==c['compound_id']};points.append((c['public_source_identifiers']['example'],obs['CDK2'],obs['CDK1']/obs['CDK2']))
fig,ax=plt.subplots(figsize=(7,4.5),layout='constrained');ax.set_xscale('log');ax.set_yscale('log')
for ex,x,y in points:
 ax.scatter(x,y,s=65,c='#007f86' if ex in [39,41] else '#777777');ax.annotate(str(ex),(x,y),xytext=(5,5 if ex!=70 else -14),textcoords='offset points',fontsize=10)
ax.axvline(100,ls='--',color='#bb5533');ax.axhline(10,ls='--',color='#bb5533');ax.set_ylim(8,240);ax.set_xlim(3,170);ax.set_xlabel('Reported CDK2/E1 IC50 (nM; lower is stronger)');ax.set_ylabel('CDK1/B1 ÷ CDK2/E1 IC50');ax.set_title('Historical analogues: potency does not settle discrimination',fontsize=12);fig.text(.5,-.02,'Patent US20240360137A1, 1 mM ATP. Labels are example numbers; no independent repeat errors supplied.',ha='center',fontsize=8);fig.savefig(out/'analogue_tradeoff.png',dpi=180,bbox_inches='tight');plt.close(fig)
designs=json.loads((R/'design-work/designs.json').read_text());im=Draw.MolsToGridImage([Chem.MolFromSmiles(d['smiles']) for d in designs],molsPerRow=2,subImgSize=(510,310),legends=[d['design_id']+' | '+('UNMEASURED hypothesis' if not d['supplied_matches'] else 'historical comparator') for d in designs]);im.save(str(out/'designs.png'))
