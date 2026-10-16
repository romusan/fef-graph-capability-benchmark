"""Corrected self-generated targets with a feasible, known 90 mm oracle.

The generator is never passed to the inverse optimizer. A failure to recover
under the stated budget remains a failure, alongside the identity sanity test.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import json,time
import numpy as np
import revision_engine as re
from revision_engine import fef,FAMILIES
from campaign import Tracker,BudgetStop,ROOT
from fast_geometry import simulate,resample,phase_alignment

def one(family,seed=901,limit=80000):
    out=ROOT/'11_validacion_completa/self_recovery_v2';out.mkdir(exist_ok=True)
    path=out/(family+'.json')
    if path.exists():return json.loads(path.read_text())
    re.configure('compact');spec=fef.PARAMS[family]
    rng=np.random.default_rng(seed);theta=np.linspace(0,2*np.pi,120,endpoint=False)
    for attempt in range(20000):
        x=re.canonical(family,rng.uniform(spec['lb'],spec['ub']))
        sim=re.simulation(family,x,theta)
        if not sim.success:continue
        width=np.ptp(sim.P[:,0])
        if width<1e-6:continue
        scale=90/width
        xx=x.copy()
        for i,name in enumerate(spec['names']):
            if name not in ['phiP','phi_c','phi_p','phase','ratio','q']:xx[i]*=scale
        if np.any(xx<spec['lb']) or np.any(xx>spec['ub']):continue
        check=re.physical_check(family,xx,1920)
        if not check['motion_valid']:continue
        sim=re.simulation(family,xx,theta);target=resample(sim.P,120)
        break
    else:raise RuntimeError(f'No feasible generator for {family}')
    truth=re.evaluate(family,xx,target,theta)
    assert truth['rms_mm']<1e-8 and truth['max_mm']<1e-8
    fef.rng=np.random.default_rng(seed+1)
    fef.N_RANDOM_PER_TOPOLOGY=900;fef.SHAPE_SEEDS_PER_TOPOLOGY=36
    fef.TOPK_RANDOM_REFINE_PER_TOPOLOGY=8;fef.RUN_DE_GLOBAL=True
    fef.DE_RESTARTS=2;fef.DE_MAXITER=55;fef.DE_POPSIZE=8
    fef.LOCAL_MAXITER=220;fef.LS_MAX_NFEV=800
    fef.ELITE_NEIGHBORHOODS=8;fef.ELITE_NEIGHBOR_SAMPLES=6
    tr=Tracker(target,[limit],re.simulation);fef.simulate_topology=tr.simulation
    try:
        while tr.count<limit:fef.search_topology(family,target,theta)
    except BudgetStop:pass
    finally:fef.simulate_topology=re.simulation
    recovered=tr.best
    recovered['physical']=re.physical_check(family,recovered['x'],1920)
    d={'family':family,'generator_seed':seed,'inverse_seed':seed+1,'generator':xx.tolist(),
       'generator_identity_rms_mm':truth['rms_mm'],'generator_identity_max_mm':truth['max_mm'],
       'generator_physical':check,'generator_given_to_optimizer':False,
       'target':target.tolist(),'simulations':tr.count,'best':recovered,'strict_found':tr.strict_best is not None}
    path.write_text(json.dumps(d,indent=2,allow_nan=False));return d

if __name__=='__main__':
    import sys
    for family in (sys.argv[1:] or FAMILIES):print(json.dumps(one(family)),flush=True)
