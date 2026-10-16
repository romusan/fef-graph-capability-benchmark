"""Reproducible descriptive diagnostics; no new optimization or seed pooling."""
import json,csv,sys
from pathlib import Path
from collections import defaultdict
import numpy as np
from summarize_v2 import ROOT,ORDER,LABELS
from render_tables_figures import latex_table
HERE=Path(__file__).parent;RESULT=HERE/'results';TABLES=ROOT/'08_manuscript/tables'
sys.path.insert(0,str(ROOT/'12_campanas_corregidas/source'))
from engineering_target import engineering_target

def main():
    targets=['easy','medium','medium_hard','bernoulli','rounded_rectangle','asymmetric_dent','naca0012_reference']
    short=['Easy','Medium','Medium-hard','Bernoulli','Rectangle','Asymmetric','NACA0012']
    descriptors=[];rows=[]
    for target,label in zip(targets,short):
        xy,meta=engineering_target(target,960)
        delta=np.roll(xy,-1,axis=0)-xy
        lengths=np.linalg.norm(delta,axis=1)
        angle=np.arctan2(delta[:,1],delta[:,0]);turn=np.angle(np.exp(1j*(np.roll(angle,-1)-angle)))
        cross=0
        for i in range(len(xy)-2):
            a=xy[i];b=xy[i+1];indices=np.arange(i+2,len(xy))
            if i==0:indices=indices[:-1]
            c=xy[indices];d=xy[(indices+1)%len(xy)]
            def area(u,v):return u[...,0]*v[...,1]-u[...,1]*v[...,0]
            cross+=int(np.sum((area(b-a,c-a)*area(b-a,d-a)<-1e-14)&(area(d-c,a-c)*area(d-c,b-c)<-1e-14)))
        desc={'target':target,'n':960,'width_mm':float(np.ptp(xy[:,0])),'height_mm':float(np.ptp(xy[:,1])),
          'polygon_length_mm':float(lengths.sum()),'absolute_turn_radians':float(np.abs(turn).sum()),
          'proper_segment_crossings':cross,'metadata':meta}
        descriptors.append(desc)
        rows.append([label,f"{desc['height_mm']:.2f}",f"{desc['polygon_length_mm']:.2f}",f"{desc['absolute_turn_radians']:.2f}"])
        np.savetxt(RESULT/f'target_{target}_960.csv',xy,delimiter=',',header='x_mm,y_mm',comments='')
    (RESULT/'target_descriptors.json').write_text(json.dumps(descriptors,indent=2))
    latex_table(TABLES/'v2_descriptors.tex','Descriptive geometry of 960-point closed target polygons at 90 mm width. Length and absolute turning are polygonal approximations. These descriptors are not a universal difficulty score; the ordered Bernoulli path additionally has a self-intersection.',
      'tab:v2_descriptors','Target & Height mm & Length mm & Absolute turn rad',rows,'@{}lrrr@{}')
    rows=[];recovery=[]
    for fam in ORDER:
        d=json.loads((ROOT/'11_validacion_completa/self_recovery_v2'/f'{fam}.json').read_text())
        b=d['best'];p=b['physical'] if b else {}
        recovery.append({'model':LABELS[fam],'identity_rms_mm':d['generator_identity_rms_mm'],
          'inverse_calls':d['simulations'],'inverse_rms_mm':b['rms_mm'] if b else None,
          'inverse_max_mm':b['max_mm'] if b else None,'any_strict':d['strict_found'],
          'retained_strict_motion':bool(b and b['strict'] and p['motion_valid']),
          'retained_quality30':bool(b and b['strict'] and p['quality30'])})
        r=recovery[-1]
        rows.append([LABELS[fam],f"{r['identity_rms_mm']:.1e}",f"{r['inverse_rms_mm']:.3f}" if b else '--',
          f"{r['inverse_max_mm']:.3f}" if b else '--',int(r['any_strict']),int(r['retained_strict_motion'])])
    (RESULT/'self_recovery_summary.json').write_text(json.dumps(recovery,indent=2))
    latex_table(TABLES/'v2_recovery.tex','Own-curve diagnostic: one feasible generator and one separate inverse search per model (80 000 simulator calls). Generator vectors are withheld from search. A denotes any encountered strict candidate; M denotes the retained best candidate being strict and motion-valid. The identity error is a consistency check, not optimization recovery.',
      'tab:v2_recovery','Model & Identity RMS & Inverse RMS & Inverse max & A & M',rows,'@{}lrrrrr@{}')
    resolution=json.loads((RESULT/'resolution_v2.json').read_text())
    assert len(resolution)==2880
    grouped=defaultdict(list)
    for d in resolution:grouped[(d['method'],d['target'],d['profile'],d['family'])].append(d)
    aggregate=[];tol=[]
    for key,dd in sorted(grouped.items()):
        assert len(dd)==30
        for n in [120,480,960]:
            rs=[next(r for r in d['resolution'] if r['n']==n) for d in dd]
            aggregate.append(dict(zip(['method','target','profile','family'],key),n=n,seeds=30,
              valid=sum(r['valid'] for r in rs),strict=sum(r['strict'] for r in rs),
              strict_motion=sum(r['strict'] and d['physical1920']['motion_valid'] for r,d in zip(rs,dd))))
        for i in range(16):
            t=dd[0]['tolerances_selected_design'][i]
            tol.append(dict(zip(['method','target','profile','family'],key),rms_threshold=t['rms_threshold'],
              max_threshold=t['max_threshold'],seeds=30,pass_count=sum(d['tolerances_selected_design'][i]['pass'] for d in dd)))
    for name,data in [('resolution_aggregate',aggregate),('tolerance_aggregate',tol)]:
        (RESULT/f'{name}.json').write_text(json.dumps(data,indent=2))
        with (RESULT/f'{name}.csv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    rows=[]
    # One descriptive row per path for the common hybrid/compact protocol: 8x30 selected vectors.
    for target,label in zip(targets,short):
        rr=[r for r in aggregate if r['method']=='hybrid' and r['profile']=='compact' and r['target']==target]
        counts=[sum(r['strict_motion'] for r in rr if r['n']==n) for n in [120,480,960]]
        rows.append([label,240,*counts])
    latex_table(TABLES/'v2_resolution.tex','Selected-vector resolution sensitivity for hybrid/compact at 20 000 calls. Counts satisfy both accuracy thresholds at each resolution and sampled motion validation at 1920 intervals. Each row describes eight distinct 30-seed model cells (240 selected vectors); it is not a pooled estimate of a model recovery probability. Invalid designs retain their denominator.',
      'tab:v2_resolution','Target & Selected vectors & N=120 & N=480 & N=960',rows,'@{}lrrrr@{}')
    print(json.dumps({'recovery':recovery,'resolution_table':rows},indent=2))

if __name__=='__main__':main()
