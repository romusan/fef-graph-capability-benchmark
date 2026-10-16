"""Generate manuscript tables/plots from complete numerical analysis only."""
import os
os.environ['MPLBACKEND']='Agg'
import json,sys,math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from summarize_v2 import LABELS,ORDER,ROOT

def latex_table(path,caption,label,header,rows,cols):
    text='\\begin{table}[H]\n\\caption{'+caption+'}\n\\label{'+label+'}\n\\centering\n\\small\n'
    text+='\\begin{tabular}{'+cols+'}\n\\toprule\n'+header+' \\\\\n\\midrule\n'
    text+='\n'.join(' & '.join(map(str,row))+' \\\\' for row in rows)
    text+='\n\\bottomrule\n\\end{tabular}\n\\end{table}\n'
    path.write_text(text,encoding='utf8')

def main():
    result=Path(__file__).parent/'results';a=json.loads((result/'analysis_v2.json').read_text())
    if not a['complete']:raise RuntimeError('Refuse final tables from incomplete data')
    paper=ROOT/'08_manuscript';tables=paper/'tables';fig=paper/'figures'
    idx={(r['method'],r['target'],r['profile'],r['family'],r['budget']):r for r in a['groups']}
    for target in ['medium','medium_hard']:
        rr=[]
        for fam in ORDER:
            c=idx[('hybrid',target,'compact',fam,20000)];r=idx[('hybrid',target,'relaxed',fam,20000)]
            rr.append([LABELS[fam],c['strict']['k'],c['strict_motion']['k'],c['strict_quality30']['k'],
              r['strict']['k'],r['strict_motion']['k'],r['strict_quality30']['k'],f"{c['median_rms']:.3f}",f"{r['median_rms']:.3f}"])
        latex_table(tables/f'v2_{target}.tex',
          f'{target.replace("_","-")} target at 20 000 simulator calls, 30 seeds per model/profile. '
          'A: any evaluated accuracy-strict candidate; M: retained strict candidate also passes sampled full-cycle motion checks; '
          'Q: that candidate additionally keeps all relevant transmission angles at least 30 degrees. '
          'RMS medians refer to minimum-RMS selected designs. C/R denote compact/relaxed. Full Wilson intervals are supplied in the machine-readable supplement.',
          'tab:v2_'+target.replace('_',''),
          'Model & A C & M C & Q C & A R & M R & Q R & RMS C & RMS R',rr,'@{}lrrrrrrrr@{}')
    rows=[]
    for fam in ORDER:
        rr=[LABELS[fam]]
        for t in ['easy','bernoulli','rounded_rectangle','asymmetric_dent']:
            r=idx[('hybrid',t,'compact',fam,20000)]
            rr += [f"{r['strict_motion']['k']}/30",f"{r['median_rms']:.3f}"]
        rows.append(rr)
    latex_table(tables/'v2_targets.tex','Additional target results at 20 000 calls and compact bounds. E: easy; B: Bernoulli; R: rounded rectangle; D: asymmetric dent. Each target/model has 30 seeds. Entries are motion-validated retained strict counts and median RMS errors in mm. Bernoulli is a published mathematical reference; the other paths are analytic, not measured industrial data.',
       'tab:v2_targets','Model & E M & E RMS & B M & B RMS & R M & R RMS & D M & D RMS',rows,'@{}lrrrrrrrr@{}')
    rows=[]
    for fam in ORDER:
        r=idx[('hybrid','naca0012_reference','compact',fam,20000)];ci=r['strict_motion']['wilson95']
        rows.append([LABELS[fam],r['strict']['k'],r['strict_motion']['k'],r['strict_quality30']['k'],f'[{ci[0]:.2f}, {ci[1]:.2f}]',f"{r['median_rms']:.3f}"])
    latex_table(tables/'v2_engineering.tex','Published NACA0012 engineering-reference contour, compact bounds, 20 000 calls and 30 seeds per model. A: accuracy strict; M: retained strict motion-valid; Q: adds all-dyad transmission quality. Wilson 95 percent intervals apply to M. Source coordinates are prescribed reference geometry, not measured mechanism trajectories.',
       'tab:v2_engineering','Model & A & M & Q & CI M & Median RMS mm',rows,'@{}lrrrrr@{}')
    rows=[]
    for fam in ORDER:
        row=[LABELS[fam]]
        for t in ['medium','medium_hard']:
            h=idx[('hybrid',t,'compact',fam,20000)];c=idx[('cma',t,'compact',fam,20000)]
            row += [h['strict_motion']['k'],c['strict_motion']['k'],f"{h['median_rms']:.3f}",f"{c['median_rms']:.3f}"]
        rows.append(row)
    latex_table(tables/'v2_independent.tex','Independent optimizer comparison under compact bounds and exactly 20 000 simulator calls per model/seed. The first four numerical columns refer to medium, the last four to medium-hard (MH). H denotes the hybrid engine and C denotes official pycma. M counts retained strict designs passing full-cycle checks; each count is out of 30. RMS is in mm.',
      'tab:v2_independent',r'& \multicolumn{4}{c}{Medium} & \multicolumn{4}{c}{Medium-hard} \\ \cmidrule(lr){2-5}\cmidrule(lr){6-9} Model & M H & M C & RMS H & RMS C & M H & M C & RMS H & RMS C',rows,'@{}lrrrrrrrr@{}')
    rows=[]
    for r in a['qd_groups']:
        if r['grid']!=24:continue
        k=r['any_strict_motion']['k'];ci=r['any_strict_motion']['wilson95']
        rows.append([r['arm'].replace('_',r'\_'),r['n'],k,f'[{ci[0]:.2f}, {ci[1]:.2f}]',r['any_quality30']['k'],f"{r['median_coverage']:.3f}",f"{r['median_qd_score']:.2f}"])
    latex_table(tables/'v2_qd.tex','Corrected QD comparison at exactly 40 000 simulator calls per arm/seed. M counts seeds with at least one retained strict motion-valid elite; Q adds the transmission filter. Coverage and QD-score use the same reference grid. Intervals are marginal Wilson 95 percent intervals.',
      'tab:v2_qd','Arm & n & M & CI & Q & Coverage & QD-score',rows,'@{}lrrrrrr@{}')
    rows=[]
    for fam in ORDER:
        p=idx[('hybrid','medium','compact',fam,20000)];n=idx[('hybrid_nosing','medium','compact',fam,20000)]
        rows.append([LABELS[fam],p['strict_motion']['k'],p['strict_quality30']['k'],n['strict_motion']['k'],n['strict_quality30']['k'],f"{p['median_rms']:.3f}",f"{n['median_rms']:.3f}"])
    latex_table(tables/'v2_penalty.tex','Controlled singularity-penalty ablation on medium/compact, 20 000 simulator calls and 30 seeds per model. The penalized row is the compatible cumulative checkpoint of the factorial, not a new independent sample. M: retained strict motion-valid; Q: additionally all transmission angles at least 30 degrees. RMS medians in mm. Removal affects both scalar cost and the LS penalty residual.',
      'tab:v2_penalty','Model & M pen & Q pen & M no & Q no & RMS pen & RMS no',rows,'@{}lrrrrrr@{}')
    ext=json.loads((result/'extended_analysis_v2.json').read_text())
    rows=[]
    for r in ext['repeated_batches']:
        ci=r['wilson95'];rows.append([r['batch'],r['n'],r['strict_motion'],f'[{ci[0]:.2f}, {ci[1]:.2f}]',f"{r['median_rms']:.3f}"])
    latex_table(tables/'v2_repetition.tex','Independent batches of CD-G/medium/compact under the same corrected 20 000-call protocol. Batch 0 uses seeds 1--30; batches 1 and 2 use 1001--1030 and 2001--2030. Main-cell counts are not inflated by pooling these diagnostic repeats. M: retained strict motion-valid; intervals are Wilson 95 percent.',
      'tab:v2_repetition','Batch & n & M & CI & Median RMS mm',rows,'@{}lrrrr@{}')
    rows=[]
    for r in a['qd_groups']:
        if r['arm']!='map1_sur1':continue
        ci=r['any_strict_motion']['wilson95'];rows.append([r['grid'],r['n'],r['any_strict_motion']['k'],f'[{ci[0]:.2f}, {ci[1]:.2f}]',f"{r['median_coverage']:.3f}",f"{r['median_qd_score']:.2f}"])
    latex_table(tables/'v2_grid.tex','Native-grid sensitivity of the full QD arm, 30 seeds each and exactly 40 000 simulator calls per run. M counts seeds with retained strict motion-valid elites. Coverage and QD-score use one common 24-by-24 reference grid for all native resolutions.',
      'tab:v2_grid','Grid side & n & M & CI & Coverage ref & QD-score ref',rows,'@{}lrrrrr@{}')
    # Budget trajectories: cumulative checkpoints, never independent samples.
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False})
    f,axes=plt.subplots(2,2,figsize=(9,6),sharex=True,sharey=True)
    colors=plt.get_cmap('tab10').colors
    for ax,(target,profile) in zip(axes.flat,[(t,p) for t in ['medium','medium_hard'] for p in ['compact','relaxed']]):
        for j,fam in enumerate(ORDER):
            rr=[idx[('hybrid',target,profile,fam,b)] for b in [10000,20000,40000]]
            ax.plot([10000,20000,40000],[r['strict_motion']['k']/30 for r in rr],'-o',ms=3,label=LABELS[fam],color=colors[j])
        ax.set_title(target.replace('_','-')+' / '+profile);ax.grid(alpha=.2);ax.set_ylim(-.02,1.02)
    for ax in axes[-1]:ax.set_xlabel('Cumulative simulator calls')
    for ax in axes[:,0]:ax.set_ylabel('Retained strict motion-valid rate')
    handles,labels=axes[0,0].get_legend_handles_labels();f.legend(handles,labels,ncol=8,loc='lower center',frameon=False)
    f.tight_layout(rect=[0,.06,1,1]);f.savefig(fig/'v2_budget.png',dpi=220);plt.close(f)
    sys.path.insert(0,str(ROOT/'12_campanas_corregidas/source'))
    from engineering_target import engineering_target as target_points
    f,axes=plt.subplots(2,4,figsize=(10,5))
    for ax,t in zip(axes.flat,['easy','medium','medium_hard','bernoulli','rounded_rectangle','asymmetric_dent','naca0012_reference']):
        xy,meta=target_points(t,960);ax.plot(xy[:,0],xy[:,1],color='#244f73');ax.set_aspect('equal');ax.set_title(t.replace('_','-'));ax.set_xlabel('mm');ax.set_ylabel('mm');ax.grid(alpha=.15)
    axes.flat[-1].set_visible(False)
    f.tight_layout();f.savefig(fig/'v2_targets.png',dpi=220);plt.close(f)
    # Export numbers as macros rather than hand-copying into the abstract.
    c=idx[('hybrid','medium_hard','compact','wattI',20000)]
    r=idx[('hybrid','medium_hard','relaxed','wattI',20000)]
    macros={'VTwoRuns':str(a['main_runs']),'VTwoQDRuns':str(a['qd_runs']),
       'CDGHardCompact':str(c['strict_motion']['k']),'CDGHardRelaxed':str(r['strict_motion']['k']),
       'CDGHardQuality':str(c['strict_quality30']['k'])}
    (tables/'v2_numbers.tex').write_text('\n'.join('\\newcommand{\\'+k+'}{'+v+'}' for k,v in macros.items()))
    model_rows=[
      ['F4','fourbar','6','One loop; tracer fixed to AB'],
      ['SC','slidercrank','5 stored; 4 shape','Crank--slider; tracer on AB'],
      ['CD-G','wattI','12','C fixed to AB; second dyad C-D-G'],
      ['RC-G','wattII','12','C fixed to rocker; second dyad C-D-G'],
      ['A-O4','stephensonI','8','Second dyad A-D-O4'],
      ['B-O2','stephensonII','8','Second dyad B-D-O2'],
      ['C-O4','stephensonIII','10','C fixed to AB; second dyad C-D-O4'],
      ['G5','gearedFiveBar','8 real + ratio','Two driven cranks; closure at C']]
    latex_table(tables/'v2_models.tex','Explicit kinematic models. Legacy identifiers are software labels, not a claim of eight distinct classical families. Tracer definitions and all bounds are provided in Methods and the appendix.',
       'tab:v2_models','Model & Legacy ID & Parameters & Construction',model_rows,'@{}lllp{5.1cm}@{}')
    sys.path.insert(0,str(ROOT/'12_campanas_corregidas/source'))
    import revision_engine as re
    import copy
    comp=copy.deepcopy(re.configure('compact'));rel=copy.deepcopy(re.configure('relaxed'))
    text=r'\small'+'\n'+r'\begin{longtable}{@{}llrrrr@{}}'+'\n'
    text+=r'\caption{Complete effective bounds in the revised experiments.}\label{tab:v2_bounds}\\'+'\n'
    text+=r'\toprule Model & Variable & Compact lower & Upper & Relaxed lower & Upper \\ \midrule\endfirsthead'+'\n'
    text+=r'\toprule Model & Variable & Compact lower & Upper & Relaxed lower & Upper \\ \midrule\endhead'+'\n'
    for fam in ORDER:
        for i,name in enumerate(comp[fam]['names']):
            values=[comp[fam]['lb'][i],comp[fam]['ub'][i],rel[fam]['lb'][i],rel[fam]['ub'][i]]
            def fmt(v):
                if abs(abs(v)-math.pi)<1e-10:return r'$-\pi$' if v<0 else r'$\pi$'
                return f'{v:.6g}'
            text+=' & '.join([LABELS[fam],r'\texttt{'+name.replace('_',r'\_')+'}']+[fmt(v) for v in values])+r' \\'+'\n'
        text+=r'\addlinespace'+'\n'
    text+=r'\bottomrule\end{longtable}\normalsize'+'\n'
    (tables/'v2_bounds.tex').write_text(text)
    print('Final tables and figures generated from complete analysis')

if __name__=='__main__':main()
