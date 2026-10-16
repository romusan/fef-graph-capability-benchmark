"""Independent trajectory audit of retained minima; no search or core edits."""
import os
os.environ['MPLBACKEND']='Agg'
import sys,json,csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
ROOT=Path(__file__).resolve().parents[1];OUT=Path(__file__).parent
sys.path.insert(0,str(ROOT/'12_campanas_corregidas/source'))
import revision_engine as re
from fast_geometry import simulate,resample,phase_alignment
from engineering_target import engineering_target as target_points

def main():
    records=[];examples={}
    for target,method in [('naca0012_reference','hybrid'),('easy','hybrid'),('medium','cma'),('medium_hard','cma')]:
        folder='engineering_results' if target.startswith('naca') else 'results'
        for fam in ['wattII','stephensonII']:
            for p in sorted((ROOT/'12_campanas_corregidas'/folder).glob(f'{target}_compact_{method}_{fam}_seed*_b20000.json')):
                d=json.loads(p.read_text());cp=next(c for c in d['checkpoints'] if c['simulations']==20000)
                b=cp['best_rms'];x=np.array(b['x']);theta=np.linspace(0,2*np.pi,120,endpoint=False)
                ok,j,margin=simulate(re.FAMILIES.index(fam),x,theta,True)
                assert ok,'Selected valid error must have a real closed trajectory'
                raw=j[:,4];y,_=target_points(target,120);pp=resample(raw,120)
                rms,mx,shift,aligned=phase_alignment(pp,y,4)
                assert abs(rms-b['rms_mm'])<1e-8
                dense_ok,dense_j,_=simulate(re.FAMILIES.index(fam),x,np.linspace(0,2*np.pi,1920,endpoint=False),True)
                singular=np.linalg.svd(raw-raw.mean(axis=0),compute_uv=False)
                ext=np.ptp(raw,axis=0);lens=np.linalg.norm(np.roll(raw,-1,axis=0)-raw,axis=1)
                physical=re.physical_check(fam,x,1920)
                r={'target':target,'method':method,'family':fam,'seed':d['seed'],'stored_rms':b['rms_mm'],
                  'recomputed_rms120':rms,'recomputed_max120':mx,'solver_valid120':bool(ok),'solver_valid1920':bool(dense_ok),'motion_valid':physical['motion_valid'],
                  'point_span_mm':ext.tolist(),'polygon_length_mm':float(lens.sum()),
                  'centered_radius_rms_mm':float(np.sqrt(np.mean(np.sum((raw-raw.mean(axis=0))**2,axis=1)))),
                  'axis_ratio':float(singular[-1]/singular[0]) if singular[0]>1e-12 else 0.,
                  'min_joint_step_mm':float(lens.min()),'max_joint_step_mm':float(lens.max()),'x':x.tolist(),
                  'physical':physical,'source':str(p.relative_to(ROOT))}
                records.append(r)
                if d['seed']==1:examples[(target,method,fam)]=(aligned,np.roll(y,shift,axis=0),raw)
    summaries=[]
    for target,method,fam in sorted(set((r['target'],r['method'],r['family']) for r in records)):
        rr=[r for r in records if (r['target'],r['method'],r['family'])==(target,method,fam)]
        assert len(rr)==30
        s={'target':target,'method':method,'family':fam,'n':30,'solver_valid':sum(r['solver_valid120'] for r in rr),
          'motion_valid':sum(r['motion_valid'] for r in rr)}
        s['solver_valid120']=sum(r['solver_valid120'] for r in rr);s['solver_valid1920']=sum(r['solver_valid1920'] for r in rr)
        for key in ['stored_rms','recomputed_rms120','polygon_length_mm','centered_radius_rms_mm','axis_ratio']:
            s[key+'_median']=float(np.median([r[key] for r in rr]))
            s[key+'_range']=[min(r[key] for r in rr),max(r[key] for r in rr)]
        s['near_line_axis_ratio_lt_0p1']=sum(r['axis_ratio']<.1 for r in rr)
        representative=min(rr,key=lambda r:abs(r['stored_rms']-s['stored_rms_median']))
        s['representative_seed']=representative['seed'];s['representative_rms']=representative['stored_rms']
        s['representative_x']=representative['x']
        x=np.array(representative['x']);ok,j,_=simulate(re.FAMILIES.index(fam),x,np.linspace(0,2*np.pi,120,endpoint=False),True)
        y,_=target_points(target,120);rms,mx,shift,aligned=phase_alignment(resample(j[:,4],120),y,4)
        examples[(target,method,fam)]=(aligned,np.roll(y,shift,axis=0),representative['seed'],rms)
        np.savetxt(OUT/f'example_{target}_{fam}_seed{representative["seed"]:03d}.csv',np.column_stack([aligned,np.roll(y,shift,axis=0)]),delimiter=',',header='path_x_mm,path_y_mm,target_x_mm,target_y_mm',comments='')
        summaries.append(s)
    ext=json.loads((ROOT/'13_cierre_editorial/results/extended_analysis_v2.json').read_text())
    counts=np.array([[r['strict_motion'],30-r['strict_motion']] for r in ext['repeated_batches']])
    chi,p,df,expected=stats.chi2_contingency(counts,correction=False)
    repetition={'table':counts.tolist(),'test':'Pearson 3x2 homogeneity; independent batches; exploratory diagnostic',
      'chi_square':float(chi),'df':int(df),'p':float(p),'minimum_expected':float(expected.min()),
      'interpretation':'Nonsignificance is not evidence of equivalence or zero batch variance.'}
    result={'records':records,'summaries':summaries,'batch_heterogeneity':repetition,
      'invalid_cost_handling':'evaluate_sim returns valid=false, cost=1e9 on closure failure, with no RMS field. Audited records have real valid trajectories.'}
    (OUT/'model_diagnostics.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    with (OUT/'model_diagnostics_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summaries[0]));w.writeheader();w.writerows(summaries)
    fig,axes=plt.subplots(2,4,figsize=(12,5.5))
    for ax,((target,method,fam),(aligned,y,seed,rms)) in zip(axes.flat,examples.items()):
        ax.plot(y[:,0],y[:,1],'k--',lw=1,label='target');ax.plot(aligned[:,0],aligned[:,1],lw=1,label='selected path')
        label={'wattII':'RC-G','stephensonII':'B-O2'}[fam]
        ax.set_title(f'{label} / {target}\n{method}, seed {seed}, RMS {rms:.3f} mm',fontsize=8);ax.set_aspect('equal');ax.grid(alpha=.2)
        ax.set_xlabel('mm');ax.set_ylabel('mm')
    axes.flat[0].legend(fontsize=7);fig.tight_layout();fig.savefig(OUT/'selected_paths_diagnosis.png',dpi=200);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(8,5.5))
    selected=[(k,v) for k,v in examples.items() if k[0] in ['naca0012_reference','medium']]
    for ax,((target,method,fam),(aligned,y,seed,rms)) in zip(axes.flat,selected):
        label={'wattII':'RC-G','stephensonII':'B-O2'}[fam]
        title='NACA0012' if target=='naca0012_reference' else 'medium'
        ax.plot(y[:,0],y[:,1],'k--',lw=1,label='Target');ax.plot(aligned[:,0],aligned[:,1],lw=1.2,label='Selected path')
        ax.set_title(f'{label} / {title}\n{method}, seed {seed}, RMS {rms:.3f} mm',fontsize=10)
        ax.set_aspect('equal');ax.grid(alpha=.2);ax.set_xlabel('mm',fontsize=9);ax.set_ylabel('mm',fontsize=9);ax.tick_params(labelsize=9)
    fig.legend(*axes.flat[0].get_legend_handles_labels(),loc='lower center',ncol=2,fontsize=9,frameon=False)
    fig.tight_layout(rect=[0,.04,1,1]);fig.savefig(OUT/'selected_paths_diagnosis_main.png',dpi=220);plt.close(fig)
    print(json.dumps({'summaries':summaries,'batch_heterogeneity':repetition},indent=2))

if __name__=='__main__':main()
