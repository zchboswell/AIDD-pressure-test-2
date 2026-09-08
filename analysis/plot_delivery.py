"""Scientific decision/evidence plots from saved tables; no hidden scoring."""
import argparse,json,pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from rdkit import Chem
from rdkit.Chem import Draw
from curate import readcsv
ROOT=pathlib.Path(__file__).resolve().parents[1]
def run(out):
 figures=out/'delivery/figures';figures.mkdir(parents=True,exist_ok=True)
 plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.titleweight':'bold','figure.facecolor':'white'})
 rows=readcsv(out/'analysis/cycle_accounting.csv');rows.sort(key=lambda r:(r['status']!='advance',r['compound_id']));colors=['#187e88' if r['status']=='advance' else '#b65a20' for r in rows];ys=np.arange(len(rows))
 fig,axes=plt.subplots(1,4,figsize=(13,5.6),sharey=True,gridspec_kw={'wspace':.22})
 configs=[('cell_TE_nM','Cellular CDK2 engagement\nIC50 (nM); lower preferred',True,150),('CDK1_nM','Recombinant CDK1/cyclin B\nIC50 (nM); context only',True,None),('solubility_uM','Kinetic solubility\nµM; higher preferred',True,None),('microsome_t_half_min','Human microsome half-life\nmin; higher preferred',False,20)]
 for ax,(key,title,log,cut) in zip(axes,configs):
  for i,r in enumerate(rows):
   value=r[key]
   if value:
    x=float(value);ax.scatter(x,i,color=colors[i],s=70,zorder=3);ax.annotate(f'{x:g}',(x,i),xytext=(5,-3),textcoords='offset points',fontsize=9)
   else:ax.text(.05,i,'unmeasured',transform=ax.get_yaxis_transform(),color='#777777',va='center',fontsize=9)
  ax.set_title(title,fontsize=10,pad=14);ax.grid(axis='x',alpha=.18);ax.set_ylim(len(rows)-.35,-.65)
  if log:
   ax.set_xscale('log');ax.xaxis.set_minor_formatter(NullFormatter())
  if key=='cell_TE_nM':
   ax.set_xlim(25,850);ax.set_xticks([50,150,500],['50','150','500'])
  elif key=='CDK1_nM':ax.set_xlim(45,5000)
  elif key=='solubility_uM':ax.set_xlim(5,250)
  else:ax.set_xlim(0,50)
  if cut:ax.axvline(cut,color='#ad253b',ls='--',lw=1,label='Project preference')
 axes[0].set_yticks(ys,[r['compound_id']+'  '+('advance' if r['status']=='advance' else 'diagnostic') for r in rows])
 fig.suptitle('Select experiments that resolve the cellular activity / CDK1 question',x=.07,ha='left',fontsize=15,fontweight='bold');fig.subplots_adjust(left=.17,right=.98,top=.79,bottom=.19)
 fig.text(.07,.065,'Supplied simulated measurements; no error bars because replicate uncertainty is unavailable. Dashed lines: program preferences.\nCDK1 and cellular engagement are different assay contexts: their ratio is not biochemical selectivity. Advances are conditional testing decisions.',fontsize=10,color='#444444')
 fig.savefig(figures/'decision_comparison.png',dpi=180);plt.close(fig)
 summary=json.loads((out/'analysis/model_summary.json').read_text());cur=json.loads((out/'analysis/curation_summary.json').read_text());oof=readcsv(out/'analysis/model_oof.csv');target=np.array([float(r['observed_pIC50']) for r in oof]);pred=np.array([float(r['rf_descriptors']) for r in oof]);legacy=json.loads(summary['original_vendor_stdout']);metrics=summary['metrics_chemical_weighted']
 fig,axes=plt.subplots(1,3,figsize=(14,4.8),gridspec_kw={'width_ratios':[1,1.3,1.2]})
 labels=list(cur['statuses']);axes[0].barh(labels,[cur['statuses'][k] for k in labels],color=['#187e88','#98a5b6','#c28e4d','#ad253b']);axes[0].invert_yaxis();axes[0].set_title('Every activity row accounted');axes[0].set_xlabel('Rows (266 total)')
 for i,k in enumerate(labels):axes[0].text(cur['statuses'][k]+2,i,str(cur['statuses'][k]),va='center')
 keys=['mean','nearest_neighbor','rf_descriptors','rf_retrospective_index'];names=['Training mean','Nearest chemical','RF: descriptors','RF: retrospective index*'];vals=[metrics[k]['mae'] for k in keys]
 axes[1].barh(names,vals,color=['#8796a5','#588a80','#187e88','#ad253b']);axes[1].invert_yaxis();axes[1].set_xlim(0,.72);axes[1].set_xlabel('Out-of-fold MAE (pIC50 log units)');axes[1].set_title('Same curated scaffold folds')
 for i,v in enumerate(vals):axes[1].text(v+.01,i,f'{v:.3f}',va='center')
 axes[2].scatter(target,pred,c='#187e88',alpha=.65,s=24);axes[2].plot([5,8.5],[5,8.5],color='#999999',ls='--');axes[2].set(xlim=(5,8.5),ylim=(5,8.5),xlabel='Historical observed pIC50',ylabel='Out-of-fold descriptor RF');axes[2].set_title(f"Descriptor RF R² = {metrics['rf_descriptors']['r2']:.3f}")
 fig.suptitle('The excellent vendor score does not establish prospective utility',x=.045,ha='left',fontsize=15,fontweight='bold');fig.subplots_adjust(left=.095,right=.98,top=.79,bottom=.28,wspace=.72)
 fig.text(.045,.075,f"Original vendor: R² {legacy['r2']:.3f}, random split of {legacy['rows']} mixed/legacy-normalized rows (a different cohort).\nCurated audit: 122 rows / 117 chemical graphs / 26 scaffold groups; folds 66, 22, 12, 11, 11. Each chemical has total weight 1.\n*Index was added after results and is absent prospectively. All models target historical biochemical export, never cellular engagement.",fontsize=10,color='#444444')
 fig.savefig(figures/'model_diagnostic.png',dpi=180);plt.close(fig)
 inventory={r['compound_id']:r for r in readcsv(out/'analysis/inventory_identity.csv') if r['identity_status']=='valid'};mols=[Chem.MolFromSmiles(inventory[r['compound_id']]['canonical_isomeric_smiles']) for r in rows]
 legends=[r['compound_id']+' | '+r['status'] for r in rows]
 image=Draw.MolsToGridImage(mols,molsPerRow=2,subImgSize=(650,350),legends=legends,useSVG=False);image.save(str(figures/'selected_chemical_structures.png'))
 receipt={'source_tables':['analysis/cycle_accounting.csv','analysis/model_oof.csv','analysis/model_summary.json','analysis/curation_summary.json','analysis/inventory_identity.csv'],'figures':['delivery/figures/decision_comparison.png','delivery/figures/model_diagnostic.png','delivery/figures/selected_chemical_structures.png'],'chemical_drawing_policy':'Full inventory isomeric graph, as supplied; no tautomer/protonation changes. CED23 stock tautomer differs from deposited ligand; no implied pose.'}
 (out/'analysis/plot_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print('Decision, model and chemical figures rendered.')
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--out',type=pathlib.Path,default=ROOT);run(a.parse_args().out.resolve())
