"""Exact simulation-budget experiments, resumable and atomic.

The budget checkpoints are cumulative trajectories, not restarted searches.
Every simulator invocation, including LS finite differences and invalid states,
is counted. Archive/effect calculations never add uncounted search evaluations.
"""
import os
for name in list(os.environ):
    if name.startswith('FEF_'):os.environ.pop(name)
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1';os.environ['NUMBA_NUM_THREADS']='1'
import sys,json,time,hashlib,contextlib,io,platform
from pathlib import Path
import numpy as np
import scipy
from scipy.optimize import differential_evolution
import revision_engine as re
from revision_engine import fef,FAMILIES
from fast_geometry import resample,phase_alignment
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

class BudgetStop(BaseException):pass

def target_points(mode,n=120):
    t=np.linspace(0,2*np.pi,4096,endpoint=False)
    if mode in ['easy','medium','medium_hard']:
        # Freeze the historically used 120-point targets for existing modes.
        return fef.generate_independent_benchmark_target(np.linspace(0,2*np.pi,n,endpoint=False),mode)
    if mode=='bernoulli':
        raw=np.column_stack([np.cos(t)/(1+np.sin(t)**2),np.sin(t)*np.cos(t)/(1+np.sin(t)**2)])
        formula='x=cos(t)/(1+sin(t)^2); y=sin(t)cos(t)/(1+sin(t)^2)'
        source='https://mathshistory.st-andrews.ac.uk/Curves/Lemniscate/'
    elif mode=='rounded_rectangle':
        # Explicit superellipse exponent 4; not described as measured data.
        raw=np.column_stack([np.sign(np.cos(t))*np.sqrt(abs(np.cos(t))),.65*np.sign(np.sin(t))*np.sqrt(abs(np.sin(t)))])
        formula='x=sign(cos(t))*sqrt(abs(cos(t))); y=.65*sign(sin(t))*sqrt(abs(sin(t)))'
        source='analytic superellipse, exponent 4, authors specified before evaluation'
    elif mode=='asymmetric_dent':
        raw=np.column_stack([np.cos(t)+.28*np.cos(3*t+.4),.7*np.sin(t)-.15*np.sin(2*t)])
        formula='x=cos(t)+.28*cos(3t+.4); y=.7*sin(t)-.15*sin(2t)'
        source='independent analytic Fourier curve, authors specified before evaluation'
    else:raise ValueError(mode)
    target=fef.normalize_target_width(resample(raw,n),90.)
    return target,{'target_mode':mode,'formula':formula,'source':source,'generated_by_family':False,'raw_samples':4096,'width_mm':90.}

def record(family,x,sim,target):
    ev=re.evaluate_sim(family,x,sim,target);sim.revision_ev=ev
    return {'family':family,'x':np.asarray(x,float).tolist(),'rms_mm':ev['rms_mm'],
            'max_mm':ev['max_mm'],'cost':ev['cost'],
            'strict':ev['accepted'],'min_delta':sim.min_delta}

class Tracker:
    def __init__(self,target,checkpoints,original):
        self.target=target;self.points=checkpoints;self.original=original
        self.count=0;self.valid=0;self.best=None;self.strict_best=None;self.by_cost=None;self.snapshots=[]
        self.candidates={}
    def simulation(self,fam,x,theta):
        if self.count>=max(self.points):raise BudgetStop()
        self.count+=1
        # Save canonical, physically meaningful vector actually simulated.
        x=re.canonical(fam,x);sim=self.original(fam,x,theta)
        if sim.success:
            self.valid+=1;r=record(fam,x,sim,self.target)
            if self.best is None or r['rms_mm']<self.best['rms_mm']:self.best=r
            if self.by_cost is None or r['cost']<self.by_cost['cost']:self.by_cost=r
            if r['strict'] and (self.strict_best is None or r['cost']<self.strict_best['cost']):self.strict_best=r
        if self.count in self.points:
            self.snapshots.append({'simulations':self.count,'valid':self.valid,'best_rms':self.best,
              'best_cost':self.by_cost,'best_strict':self.strict_best,'any_strict':self.strict_best is not None})
        return sim

def run_one(job,outdir=None):
    family,seed,profile,mode,method,limit=job
    name=f'{mode}_{profile}_{method}_{family}_seed{seed:03d}_b{limit}'
    out=Path(outdir) if outdir else ROOT/'results'
    out.mkdir(parents=True,exist_ok=True);path=out/(name+'.json')
    source_hash=hashlib.sha256(b''.join((HERE/n).read_bytes() for n in
        ['historical_engine.py','fast_geometry.py','revision_engine.py','campaign.py'])).hexdigest()
    if path.exists():
        d=json.loads(path.read_text())
        if d['source_hash']!=source_hash:raise RuntimeError(f'Stale results cannot be reused: {path}')
        return {'name':name,'cached':True,'seconds':d['seconds']}
    start=time.time();re.configure(profile)
    fef.rng=np.random.default_rng(seed);np.random.seed(seed)
    fef.N_RANDOM_PER_TOPOLOGY=900;fef.SHAPE_SEEDS_PER_TOPOLOGY=36
    fef.TOPK_RANDOM_REFINE_PER_TOPOLOGY=8;fef.RUN_DE_GLOBAL=True
    fef.DE_RESTARTS=2;fef.DE_MAXITER=55;fef.DE_POPSIZE=8
    fef.LOCAL_MAXITER=220;fef.LS_MAX_NFEV=800;fef.RUN_LS_POLISH=True
    fef.ELITE_NEIGHBORHOODS=8;fef.ELITE_NEIGHBOR_SAMPLES=6;fef.RANDOM_SAMPLER='lhs'
    fef.W_SING=0. if method=='hybrid_nosing' else .08
    target,meta=target_points(mode);theta=np.linspace(0,2*np.pi,120,endpoint=False)
    points=sorted(set([limit//4,limit//2,limit]))
    tr=Tracker(target,points,re.simulation);fef.simulate_topology=tr.simulation
    rng=np.random.default_rng(seed)
    spec=fef.PARAMS[family];lb,ub=spec['lb'],spec['ub']
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            if method in ['hybrid','hybrid_nosing']:
                while tr.count<limit:fef.search_topology(family,target,theta)
            elif method=='cma':
                import cma
                while tr.count<limit:
                    # Independent CMA starts use normalized box, no FEF templates.
                    es=cma.CMAEvolutionStrategy(rng.uniform(.1,.9,len(lb)),.25,
                        {'bounds':[0,1],'seed':seed+tr.count+1,'verbose':-9,'verb_log':0})
                    while not es.stop() and tr.count<limit:
                        xx=es.ask();costs=[]
                        for u in xx:costs.append(fef.evaluate_candidate(family,lb+np.asarray(u)*(ub-lb),target,theta)['cost'])
                        es.tell(xx,costs)
            else:raise ValueError(method)
    except BudgetStop:pass
    finally:fef.simulate_topology=re.simulation
    assert tr.count==limit
    # Only post-hoc verification: this adds no search candidates.
    for snap in tr.snapshots:
        for field in ['best_rms','best_cost','best_strict']:
            r=snap[field]
            if r is not None:
                r['physical']=re.physical_check(family,r['x'],1920)
    d={'name':name,'protocol':'revision-v2-cumulative-simulations','source_hash':source_hash,
       'family':family,'seed':seed,'profile':profile,'target':meta,'method':method,
       'limit':limit,'simulations':tr.count,'valid':tr.valid,'checkpoints':tr.snapshots,
       'bounds':{'names':spec['names'],'lb':lb.tolist(),'ub':ub.tolist()},
       'gear_ratios':re.GEAR_RATIOS.tolist() if family=='gearedFiveBar' else None,
       'versions':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
       'seconds':time.time()-start}
    tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(d,indent=2,allow_nan=False));tmp.replace(path)
    return {'name':name,'seconds':d['seconds'],'best_rms':tr.best['rms_mm'] if tr.best else None,
            'strict':tr.strict_best is not None}

if __name__=='__main__':
    print(json.dumps(run_one((sys.argv[1],int(sys.argv[2]),sys.argv[3],sys.argv[4],sys.argv[5],int(sys.argv[6])))),flush=True)
