"""Re-evaluate archived selected designs without clipping, projection or repair."""
import os,sys,json,csv,re,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'12_campanas_corregidas/source'))
import numpy as np
from fast_geometry import simulate,resample,phase_alignment
from revision_engine import FAMILIES,physical_check,fef

def main():
    data=ROOT.parents[1]/'repo_codigo_nuevo/datos'
    if not data.exists():data=ROOT.parents[1]/'datos'
    sources=[]
    for tag,folder in [('easy',data/'benchmark/easy_all_families'),
      ('medium_initial',data/'benchmark/medium_wattI'),('medium_all8',data/'policy/all_families_medium'),
      ('mh_historical_wattI_compact',data/'packaging/medium_hard_compact_wattI'),
      ('mh_historical_relaxed',data/'benchmark/medium_hard_relaxed'),
      ('mh_historical_geared',data/'packaging/medium_hard_geared_five_bar'),
      ('mh_new_compact',ROOT/'03_target_benchmark/medium_hard_compact'),
      ('mh_new_relaxed',ROOT/'03_target_benchmark/medium_hard_relaxed'),
      ('slider180_medium',ROOT/'04_packaging_singularity/slidercrank_180_medium')]:
        for p in sorted(folder.glob('*/candidate_archive.csv')):
            if tag=='medium_initial' and not p.parent.name.startswith('paper_stats_run_wattI_medium_seed'):continue
            with p.open(encoding='utf-8-sig',newline='') as f:rr=list(csv.DictReader(f))
            family_default=p.parent.name.split('_seed')[0]
            groups={}
            for r in rr:
                fam=r.get('topology',family_default)
                if fam not in FAMILIES and tag=='slider180_medium':fam='slidercrank'
                if fam not in groups or float(r['rms_mm'])<float(groups[fam]['rms_mm']):groups[fam]=r
            for fam,r in groups.items():
                names=fef.PARAMS[fam]['names']
                x=np.array([float(r[k]) for k in names])
                actual_tag=tag
                if tag=='mh_historical_geared':
                    actual_tag=tag+('_compact' if '_compact_' in p.parent.name else '_relaxed')
                sources.append((actual_tag,p,fam,x,float(r['rms_mm']),float(r['max_mm'])))
    result=[];out=Path(__file__).parent
    for i,(tag,p,fam,x,rms,mx) in enumerate(sources):
        checks=[physical_check(fam,x,n,legacy=True) for n in [120,480,1920]]
        orig_strict=rms<=1 and mx<=2
        dense={}
        if orig_strict:
            mode='easy' if tag=='easy' else 'medium_hard' if tag.startswith('mh_') else 'medium'
            theta=np.linspace(0,2*np.pi,960,endpoint=False)
            ok,j,_=simulate(FAMILIES.index(fam),x,theta,False)
            if ok:
                target,_=fef.generate_independent_benchmark_target(theta,mode)
                dr,dm,_,_=phase_alignment(resample(j[:,4],960),target,1)
                dense={'rms':dr,'max':dm,'strict':bool(dr<=1 and dm<=2)}
        rec={'campaign':tag,'path':str(p),'family':fam,'x':x.tolist(),
           'seed':int(re.search(r'seed(\d+)',p.parent.name).group(1)),
           'original_rms':rms,'original_max':mx,'original_strict':orig_strict,
           'checks':checks,'dense_accuracy':dense,
           'motion_valid_all_grids':all(c['motion_valid'] for c in checks),
           'quality30_all_grids':all(c['quality30'] for c in checks)}
        result.append(rec)
        if (i+1)%50==0:print(f'validated {i+1}/{len(sources)}',flush=True)
    qd=[]
    for p in sorted((data/'policy/qd_arms').glob('*.json')):
        d=json.loads(p.read_text())
        if d['config'] not in ['Estrict','Bplus','Bmatched']:continue
        for i,r in enumerate(d['refined']):
            if not r['strict']:continue
            x=r.get('x')
            qd.append({'path':str(p),'index':i,'family':r['topology'],'has_vector':x is not None,
                'check':physical_check(r['topology'],x,1920,True) if x is not None else None})
    (out/'archived_validation.json').write_text(json.dumps({'designs':result,'qd':qd},indent=2),encoding='utf8')
    summarize(result,out)

def summarize(result,out):
    # Exclude partially archived exploratory medium-hard directories sharing
    # the parent folder of the 30-seed medium campaign.
    result=[r for r in result if r['campaign']!='medium_initial' or
        Path(r['path']).parent.name.startswith('paper_stats_run_wattI_medium_seed')]
    for r in result:
        if r['campaign']=='mh_historical_geared':
            r['campaign']+=('_compact' if '_compact_' in Path(r['path']).parent.name else '_relaxed')
    doc=json.loads((out/'archived_validation.json').read_text(encoding='utf8'))
    doc['designs']=result
    (out/'archived_validation.json').write_text(json.dumps(doc,indent=2),encoding='utf8')
    summary={}
    for r in result:
        key=r['campaign']+'/'+r['family']
        s=summary.setdefault(key,{'n':0,'accuracy':0,'accuracy_motion':0,'accuracy_quality30':0,'accuracy_dense':0})
        s['n']+=1;s['accuracy']+=r['original_strict']
        s['accuracy_motion']+=r['original_strict'] and r['motion_valid_all_grids']
        s['accuracy_quality30']+=r['original_strict'] and r['quality30_all_grids']
        s['accuracy_dense']+=r.get('dense_accuracy',{}).get('strict',False) and r['motion_valid_all_grids']
    (out/'archived_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':
    if '--summarize' in sys.argv:
        out=Path(__file__).parent;summarize(json.loads((out/'archived_validation.json').read_text())['designs'],out)
    else:main()
