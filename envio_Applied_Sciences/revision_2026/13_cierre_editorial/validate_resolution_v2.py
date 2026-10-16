"""Resolution and tolerance sensitivity of selected corrected designs.

Post-hoc only: no optimization at the dense grids, no replacement of failed
seeds. Error tolerances apply to the minimum-RMS selected vector, explicitly
distinct from the any-candidate search indicator used for budget curves.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'12_campanas_corregidas/source'))
from campaign import FAMILIES
from engineering_target import engineering_target as target_points
from fast_geometry import simulate,resample,phase_alignment
from revision_engine import physical_check

def main():
    out=Path(__file__).parent/'results';out.mkdir(exist_ok=True)
    pp=list((ROOT/'12_campanas_corregidas/results').glob('*.json'))+list((ROOT/'12_campanas_corregidas/engineering_results').glob('*.json'))
    if len(pp)!=2880:raise RuntimeError('Formal campaigns must finish before resolution validation')
    records=[]
    for i,p in enumerate(sorted(pp)):
        d=json.loads(p.read_text());cp=next((c for c in d['checkpoints'] if c['simulations']==20000),None)
        if cp is None:raise RuntimeError('Missing base checkpoint')
        b=cp['best_rms']
        if b is None:
            records.append({'name':d['name'],'valid':False});continue
        x=np.array(b['x']);fam=d['family'];mode=d['target']['target_mode']
        vals=[]
        for n in [120,480,960]:
            theta=np.linspace(0,2*np.pi,n,endpoint=False)
            ok,j,_=simulate(FAMILIES.index(fam),x,theta,True)
            if not ok:vals.append({'n':n,'valid':False,'strict':False});continue
            target,_=target_points(mode,n)
            rms,mx,_,_=phase_alignment(resample(j[:,4],n),target,1)
            vals.append({'n':n,'valid':True,'rms':rms,'max':mx,'strict':bool(rms<=1 and mx<=2)})
        tol=[]
        for r in [.5,1.,1.5,2.]:
            for m in [1.,2.,3.,4.]:tol.append({'rms_threshold':r,'max_threshold':m,'pass':bool(vals[0].get('rms',float('inf'))<=r and vals[0].get('max',float('inf'))<=m)})
        records.append({'name':d['name'],'family':fam,'profile':d['profile'],'method':d['method'],'seed':d['seed'],
          'target':mode,'budget_selected':20000,'resolution':vals,'tolerances_selected_design':tol,
          'physical1920':physical_check(fam,x,1920)})
        if (i+1)%100==0:print(f'resolution {i+1}/{len(pp)}',flush=True)
    (out/'resolution_v2.json').write_text(json.dumps(records,indent=2,allow_nan=False))
    print('Resolution validation complete',flush=True)

if __name__=='__main__':main()
