"""Seed-blocked component, optimizer, penalty, grid and repetition summaries."""
import json,sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy import stats
from summarize_v2 import ROOT,LABELS,ORDER,holm,wilson
HERE=Path(__file__).parent;CAMP=ROOT/'12_campanas_corregidas'

def paired_binary(a,b):
    ss=sorted(set(a)&set(b))
    if len(ss)!=30:raise ValueError('A paired comparison requires exactly the 30 defined blocks')
    d=np.array([int(a[s])-int(b[s]) for s in ss],float)
    gains=int((d>0).sum());losses=int((d<0).sum())
    boot=np.random.default_rng(16016).choice(d,(20000,30),replace=True).mean(axis=1)
    return {'n':30,'difference_a_minus_b':float(d.mean()),'difference_boot95':np.percentile(boot,[2.5,97.5]).tolist(),
      'a_only':gains,'b_only':losses,
      'mcnemar_p':float(stats.binomtest(gains,gains+losses,.5).pvalue) if gains+losses else 1.}

def paired_continuous(a,b):
    ss=sorted(set(a)&set(b));d=np.array([a[s]-b[s] for s in ss])
    if len(ss)!=30:raise ValueError('Missing blocks')
    boot=np.random.default_rng(16017).choice(d,(20000,30),replace=True).mean(axis=1)
    return {'n':30,'mean_difference_a_minus_b':float(d.mean()),'mean_difference_boot95':np.percentile(boot,[2.5,97.5]).tolist(),
       'wilcoxon_p':float(stats.wilcoxon(d).pvalue) if np.any(d) else 1.}

def adjust(items,field):
    for item,adj in zip(items,holm([i[field] for i in items])):item[field+'_holm']=adj

def main():
    result=HERE/'results';analysis=json.loads((result/'analysis_v2.json').read_text())
    if not analysis['complete']:raise RuntimeError('Incomplete evaluation')
    docs=[json.loads(p.read_text()) for p in (CAMP/'results').glob('*.json')]
    idx={}
    for d in docs:
        for cp in d['checkpoints']:
            if cp['simulations']!=20000:continue
            r=cp['best_rms'];s=cp['best_strict']
            if r is None:raise ValueError('No valid selected RMS; do not replace missing errors by a finite sentinel')
            idx[(d['method'],d['target']['target_mode'],d['profile'],d['family'],d['seed'])]={
               'rms':r['rms_mm'],
               'motion':bool(s and s['physical']['motion_valid']),
               'quality':bool(s and s['physical']['quality30']),
               'angle':s['physical'].get('min_angle_deg') if s else None}
    independent=[];penalty=[]
    for target in ['medium','medium_hard']:
        for fam in ORDER:
            h={s:idx[('hybrid',target,'compact',fam,s)] for s in range(1,31)}
            c={s:idx[('cma',target,'compact',fam,s)] for s in range(1,31)}
            independent.append({'target':target,'family':fam,
              **paired_binary({s:h[s]['motion'] for s in h},{s:c[s]['motion'] for s in c}),
              'rms':paired_continuous({s:h[s]['rms'] for s in h},{s:c[s]['rms'] for s in c})})
    adjust(independent,'mcnemar_p')
    for rec,adj in zip(independent,holm([r['rms']['wilcoxon_p'] for r in independent])):rec['rms']['wilcoxon_p_holm']=adj
    for fam in ORDER:
        h={s:idx[('hybrid','medium','compact',fam,s)] for s in range(1,31)}
        c={s:idx[('hybrid_nosing','medium','compact',fam,s)] for s in range(1,31)}
        penalty.append({'family':fam,'direction':'penalized minus no penalty',
           'motion':paired_binary({s:h[s]['motion'] for s in h},{s:c[s]['motion'] for s in c}),
           'quality':paired_binary({s:h[s]['quality'] for s in h},{s:c[s]['quality'] for s in c}),
           'rms':paired_continuous({s:h[s]['rms'] for s in h},{s:c[s]['rms'] for s in c})})
    for field,key in [('motion','mcnemar_p'),('quality','mcnemar_p'),('rms','wilcoxon_p')]:
        for r,adj in zip(penalty,holm([p[field][key] for p in penalty])):r[field][key+'_holm']=adj
    qd=defaultdict(dict)
    for p in (CAMP/'qd_results').glob('*.json'):
        d=json.loads(p.read_text())
        if d['seed'] not in range(1,31):continue
        qd[(d['arm'],d['native_grid'])][d['seed']]=d
    tests=[]
    pairs=[('map effect without surrogate','map1_sur0','map0_sur0'),
           ('map effect with surrogate','map1_sur1','map0_sur1'),
           ('surrogate effect without map','map0_sur1','map0_sur0'),
           ('surrogate effect with map','map1_sur1','map1_sur0'),
           ('full QD versus capped hybrid','map1_sur1','baseline')]
    for name,armA,armB in pairs:
        a,b=qd[(armA,24)],qd[(armB,24)]
        tests.append({'contrast':name,**paired_binary({s:r['any_strict_motion'] for s,r in a.items()},{s:r['any_strict_motion'] for s,r in b.items()}),
          'coverage':paired_continuous({s:r['coverage_reference'] for s,r in a.items()},{s:r['coverage_reference'] for s,r in b.items()})})
    adjust(tests,'mcnemar_p')
    for r,adj in zip(tests,holm([r['coverage']['wilcoxon_p'] for r in tests])):r['coverage']['wilcoxon_p_holm']=adj
    grid=[]
    for side in [12,18,36]:
        a,b=qd[('map1_sur1',side)],qd[('map1_sur1',24)]
        grid.append({'native_grid':side,'reference_grid':24,
          **paired_binary({s:r['any_strict_motion'] for s,r in a.items()},{s:r['any_strict_motion'] for s,r in b.items()}),
          'coverage':paired_continuous({s:r['coverage_reference'] for s,r in a.items()},{s:r['coverage_reference'] for s,r in b.items()})})
    adjust(grid,'mcnemar_p')
    for r,adj in zip(grid,holm([r['coverage']['wilcoxon_p'] for r in grid])):r['coverage']['wilcoxon_p_holm']=adj
    repetition={0:[idx[('hybrid','medium','compact','wattI',s)] for s in range(1,31)]}
    for p in (CAMP/'repeatability_results').glob('*.json'):
        d=json.loads(p.read_text());s=next(c for c in d['checkpoints'] if c['simulations']==20000)['best_strict'];b=d['checkpoints'][-1]['best_rms']
        batch=d['seed']//1000
        repetition.setdefault(batch,[]).append({'motion':bool(s and s['physical']['motion_valid']),'rms':b['rms_mm']})
    rep=[]
    for batch,rr in sorted(repetition.items()):
        if len(rr)!=30:raise ValueError('Incomplete repetition batch')
        k=sum(r['motion'] for r in rr)
        rep.append({'batch':batch,'n':30,'strict_motion':k,'wilson95':wilson(k,30),
         'median_rms':float(np.median([r['rms'] for r in rr]))})
    result_doc={'independent_optimizer':independent,'penalty_ablation':penalty,
        'qd_component_contrasts':tests,'qd_grid_contrasts':grid,'repeated_batches':rep,
        'between_batch_rate_sample_variance':float(np.var([r['strict_motion']/30 for r in rep],ddof=1)),
        'multiplicity':'Separate Holm families for each stated contrast block and endpoint. Bootstrap intervals are marginal.',
        'seed_dependence':'Seed blocks are experimental replicates; elite vectors and budget checkpoints are not extra sample units.'}
    (result/'extended_analysis_v2.json').write_text(json.dumps(result_doc,indent=2,allow_nan=False))
    print('Extended comparisons complete')

if __name__=='__main__':main()
