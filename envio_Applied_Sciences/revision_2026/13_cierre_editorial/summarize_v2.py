"""Single source of numerical tables; refuses incomplete formal campaigns."""
import os
os.environ['MPLBACKEND']='Agg'
import sys,json,csv,math,itertools
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy import stats

ROOT=Path(__file__).resolve().parents[1];CAMP=ROOT/'12_campanas_corregidas'
LABELS={'fourbar':'F4','slidercrank':'SC','wattI':'CD-G','wattII':'RC-G',
 'stephensonI':'A-O4','stephensonII':'B-O2','stephensonIII':'C-O4','gearedFiveBar':'G5'}
ORDER=list(LABELS)

def wilson(k,n):
    if n==0:return [None,None]
    z=stats.norm.ppf(.975);p=k/n;den=1+z*z/n
    mid=(p+z*z/(2*n))/den;rad=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [float(max(0,mid-rad)),float(min(1,mid+rad))]

def holm(p):
    order=np.argsort(p);adj=np.empty(len(p));prev=0.
    for i,j in enumerate(order):prev=max(prev,min(1.,(len(p)-i)*p[j]));adj[j]=prev
    return adj.tolist()

def main():
    partial='--allow-partial' in sys.argv
    out=Path(__file__).parent/('partial' if partial else 'results');out.mkdir(exist_ok=True)
    rows=[];sources=[]
    main_files=sorted((CAMP/'results').glob('*.json'))
    engineering_files=sorted((CAMP/'engineering_results').glob('*.json'))
    for p in main_files+engineering_files:
        d=json.loads(p.read_text());sources.append(str(p.relative_to(ROOT)))
        for cp in d['checkpoints']:
            b=cp['best_rms'];s=cp['best_strict']
            rows.append({'target':d['target']['target_mode'],'profile':d['profile'],'method':d['method'],
              'family':d['family'],'model':LABELS[d['family']],'seed':d['seed'],'budget':cp['simulations'],
              'rms':b['rms_mm'] if b else None,'max':b['max_mm'] if b else None,
              'strict':cp['any_strict'],'strict_motion':bool(s and s['physical']['motion_valid']),
              'strict_quality30':bool(s and s['physical']['quality30']),
              'selected_strict':s,'best_rms':b})
    expected=960+240+480+720+240+240
    qd=[json.loads(p.read_text()) for p in sorted((CAMP/'qd_results').glob('*.json'))
         if int(p.stem.split('seed')[-1]) in range(1,31)]
    complete=len(main_files)==2640 and len(engineering_files)==240 and len(qd)==240
    if not partial and not complete:raise RuntimeError(f'Incomplete experiments: main={len(sources)}/{expected}, QD={len(qd)}/240')
    groups=defaultdict(list)
    for r in rows:groups[tuple(r[k] for k in ['method','target','profile','family','budget'])].append(r)
    summaries=[]
    for key,rr in sorted(groups.items()):
        seeds=[r['seed'] for r in rr]
        if len(set(seeds))!=len(seeds):raise ValueError(f'Duplicate seeds {key}')
        if not partial and sorted(seeds)!=list(range(1,31)):raise ValueError(f'Missing seeds {key}')
        rec=dict(zip(['method','target','profile','family','budget'],key));rec['model']=LABELS[rec['family']];rec['n']=len(rr)
        for metric in ['strict','strict_motion','strict_quality30']:
            k=sum(r[metric] for r in rr);rec[metric]={'k':k,'rate':k/len(rr),'wilson95':wilson(k,len(rr))}
        rms=[r['rms'] for r in rr if r['rms'] is not None]
        rec['median_rms']=float(np.median(rms)) if rms else None
        rec['iqr_rms']=np.percentile(rms,[25,75]).tolist() if rms else None
        summaries.append(rec)
    # Predeclared profile contrasts at B, by target/model: 16 paired RMS tests.
    contrasts=[]
    for target in ['medium','medium_hard']:
        for family in ORDER:
            c={r['seed']:r for r in rows if (r['method'],r['target'],r['profile'],r['family'],r['budget'])==('hybrid',target,'compact',family,20000)}
            a={r['seed']:r for r in rows if (r['method'],r['target'],r['profile'],r['family'],r['budget'])==('hybrid',target,'relaxed',family,20000)}
            ss=sorted(set(c)&set(a))
            if len(ss)<30:continue
            diff=np.array([c[s]['rms']-a[s]['rms'] for s in ss])
            p=float(stats.wilcoxon(diff).pvalue) if np.any(diff) else 1.
            boot=np.random.default_rng(497).choice(diff,(20000,len(diff)),replace=True).mean(axis=1)
            gains=sum(c[s]['strict_motion'] and not a[s]['strict_motion'] for s in ss)
            losses=sum(a[s]['strict_motion'] and not c[s]['strict_motion'] for s in ss)
            contrasts.append({'target':target,'family':family,'n':len(ss),'mean_rms_comp_minus_rel':float(diff.mean()),
              'mean_shift_boot95':np.percentile(boot,[2.5,97.5]).tolist(),
              'p_rms_wilcoxon':p,'compact_only_motion':gains,'relaxed_only_motion':losses,
              'p_motion_mcnemar':float(stats.binomtest(gains,gains+losses,.5).pvalue) if gains+losses else 1.})
    if contrasts:
        # Separate families for errors and pass indicators; each has 16 tests.
        for field in ['p_rms_wilcoxon','p_motion_mcnemar']:
            for r,adj in zip(contrasts,holm([r[field] for r in contrasts])):r[field+'_holm']=adj
    qgroups=defaultdict(list)
    for d in qd:qgroups[(d['arm'],d['native_grid'])].append(d)
    qsummary=[]
    for (arm,grid),dd in sorted(qgroups.items()):
        assert all(d['simulations']==40000 for d in dd)
        r={'arm':arm,'grid':grid,'n':len(dd)}
        for metric in ['any_strict','any_strict_motion','any_quality30']:
            k=sum(d[metric] for d in dd);r[metric]={'k':k,'wilson95':wilson(k,len(dd))}
        r['median_coverage']=float(np.median([d['coverage_reference'] for d in dd]))
        r['median_qd_score']=float(np.median([d['qd_score_reference'] for d in dd]));qsummary.append(r)
    result={'complete':complete,'main_runs':len(sources),'expected_main_runs':expected,'qd_runs':len(qd),
       'groups':summaries,'profile_contrasts':contrasts,'qd_groups':qsummary,'sources':sources,
       'interpretation':'Budget checkpoints are dependent cumulative observations. Motion/quality refer to retained selected designs.'}
    (out/'analysis_v2.json').write_text(json.dumps(result,indent=2,allow_nan=False))
    with (out/'per_seed_v2.csv').open('w',newline='',encoding='utf8') as f:
        fields=[k for k in rows[0] if k not in ['selected_strict','best_rms']] if rows else []
        w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
    print(json.dumps({'complete':complete,'main_runs':len(sources),'qd_runs':len(qd),'groups':len(summaries)}))

if __name__=='__main__':main()
