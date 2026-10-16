"""Post-hoc sampled point envelope of the 960 selected factorial base designs."""
import os
os.environ['MPLBACKEND']='Agg'
import sys,json,csv
from pathlib import Path
from collections import defaultdict
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=Path(__file__).parent
sys.path.insert(0,str(ROOT/'12_campanas_corregidas/source'))
from fast_geometry import simulate
from revision_engine import FAMILIES
sys.path.insert(0,str(ROOT/'13_cierre_editorial'))
from render_tables_figures import latex_table
from summarize_v2 import ORDER,LABELS

def main():
    records=[]
    paths=[p for p in (ROOT/'12_campanas_corregidas/results').glob('*_b40000.json')]
    assert len(paths)==960
    for p in sorted(paths):
        d=json.loads(p.read_text());cp=next(c for c in d['checkpoints'] if c['simulations']==20000)
        b=cp['best_rms'];x=np.array(b['x']);fam=d['family']
        ok,j,_=simulate(FAMILIES.index(fam),x,np.linspace(0,2*np.pi,1921),True)
        r={k:d[k] for k in ['family','profile','seed']};r['target']=d['target']['target_mode'];r['valid1920']=bool(ok)
        if ok:
            anchor=np.array([x[4],0.]) if fam=='slidercrank' else np.zeros(2)
            points=j.reshape(-1,2)
            anchors=np.array([anchor]) if fam=='slidercrank' else np.array([[0.,0.],[x[0],0.]])
            if fam in ['wattI','wattII']:anchors=np.vstack([anchors,[x[6],x[7]]])
            points=np.vstack([points,anchors]);lo=points.min(axis=0);hi=points.max(axis=0);wh=hi-lo
            assert np.all(wh>=0) and np.isfinite(points).all()
            r.update(width_mm=float(wh[0]),height_mm=float(wh[1]),area_mm2=float(np.prod(wh)),
              maximum_anchor_radius_mm=float(np.linalg.norm(points-anchor,axis=1).max()),min_xy=lo.tolist(),max_xy=hi.tolist())
        records.append(r)
    groups=defaultdict(list)
    for r in records:groups[(r['target'],r['family'],r['profile'])].append(r)
    summaries=[]
    for key,rr in sorted(groups.items()):
        assert len(rr)==30;valid=[r for r in rr if r['valid1920']]
        s=dict(zip(['target','family','profile'],key),n=30,valid1920=len(valid))
        for metric in ['area_mm2','maximum_anchor_radius_mm','width_mm','height_mm']:
            s[metric+'_median']=float(np.median([r[metric] for r in valid])) if valid else None
        summaries.append(s)
    (OUT/'envelope_diagnostics.json').write_text(json.dumps({'records':records,'summaries':summaries},indent=2,allow_nan=False))
    with (OUT/'envelope_diagnostics_summary.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summaries[0]));w.writeheader();w.writerows(summaries)
    idx={(s['target'],s['family'],s['profile']):s for s in summaries};rows=[]
    for target,label in [('medium','M'),('medium_hard','MH')]:
        for fam in ORDER:
            c=idx[(target,fam,'compact')];r=idx[(target,fam,'relaxed')]
            rows.append([label,LABELS[fam],c['valid1920'],r['valid1920'],f"{c['area_mm2_median']/1000:.2f}" if c['valid1920'] else '--',f"{r['area_mm2_median']/1000:.2f}" if r['valid1920'] else '--',f"{c['maximum_anchor_radius_mm_median']:.1f}" if c['valid1920'] else '--',f"{r['maximum_anchor_radius_mm_median']:.1f}" if r['valid1920'] else '--'])
    latex_table(ROOT/'08_manuscript/tables/v3_envelopes.tex',
      'Sampled point envelopes of minimum-RMS selected designs at 20 000 calls. M/MH: medium/medium-hard; C/R: compact/relaxed. Valid counts are out of 30 and refer to geometric closure on the 1920-interval grid. Area is in thousands of square millimetres; radius is in millimetres. Medians exclude failed dense closures and keep their counts visible. These descriptive point envelopes do not model body thickness or collisions.',
      'tab:v3_envelopes','Target & Model & Valid C & Valid R & Area C & Area R & Radius C & Radius R',rows,'@{}llrrrrrr@{}')
    print('Envelopes checked:',len(records))

if __name__=='__main__':main()
