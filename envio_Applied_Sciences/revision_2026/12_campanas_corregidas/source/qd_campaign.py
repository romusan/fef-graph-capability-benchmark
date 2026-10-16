"""Budget-exact QD rerun with controlled MAP-Elites and surrogate factors.

Uses the original weighted-local-ridge surrogate, descriptor definitions and
Gaussian mutation. All candidate vectors are now retained for symmetric QA.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,json,time,hashlib,contextlib,io
from pathlib import Path
import numpy as np
import revision_engine as re
from campaign import Tracker,BudgetStop,target_points,ROOT,HERE
from revision_engine import fef,FAMILIES

PROJECT=ROOT.parents[1]
sys.path.insert(0,str(PROJECT/'qd_experiments'))
sys.modules['fef_mechanism_family_synthesis']=fef
import qd_core as qc
import qd_configs as cfg

class QDTracker(Tracker):
    def __init__(self,target,limit):
        super().__init__(target,[limit],re.simulation)
        self.reference={};self.strict_archive={};self.pool={};self.cap=limit
    def simulation(self,fam,x,theta):
        if self.count>=self.cap:raise BudgetStop()
        sim=super().simulation(fam,x,theta)
        if sim.success:
            ev=sim.revision_ev
            desc=qc.descriptors(ev);ev['_d']=desc
            key=qc.cell_of(*desc)
            r={'topology':fam,'x':ev['x'].copy(),'rms_mm':ev['rms_mm'],'max_mm':ev['max_mm'],'_d':desc}
            inc=self.reference.get(key)
            if inc is None or r['rms_mm']<inc['rms_mm']:self.reference[key]=r
            if ev['accepted']:
                # Strict best per family/reference cell, separate from RMS elite.
                sk=(fam,*key);inc=self.strict_archive.get(sk)
                if inc is None or r['rms_mm']<inc['rms_mm']:self.strict_archive[sk]=r
        return sim

def run_qd(job):
    seed,map_elites,surrogate,grid,baseline=job
    label='baseline' if baseline else f'map{int(map_elites)}_sur{int(surrogate)}'
    out=ROOT/'12_campanas_corregidas/qd_results';out.mkdir(exist_ok=True)
    path=out/f'{label}_grid{grid}_seed{seed:03d}.json'
    sources=[HERE/n for n in ['historical_engine.py','fast_geometry.py','revision_engine.py','campaign.py','qd_campaign.py']]
    sources += [PROJECT/'qd_experiments'/n for n in ['qd_core.py','qd_configs.py']]
    sha=hashlib.sha256(b''.join(p.read_bytes() for p in sources)).hexdigest()
    if path.exists():
        d=json.loads(path.read_text())
        if d['source_hash']!=sha:raise RuntimeError('QD source mismatch')
        return {'name':path.name,'cached':True}
    start=time.time();re.configure('compact');fef.rng=np.random.default_rng(seed+7)
    fef.W_SING=.08;fef.N_RANDOM_PER_TOPOLOGY=120;fef.SHAPE_SEEDS_PER_TOPOLOGY=16
    fef.RUN_DE_GLOBAL=False;fef.TOPK_RANDOM_REFINE_PER_TOPOLOGY=4
    fef.LOCAL_MAXITER=100;fef.LS_MAX_NFEV=300;fef.RUN_LS_POLISH=True
    fef.ELITE_NEIGHBORHOODS=6;fef.ELITE_NEIGHBOR_SAMPLES=4
    qc.GRID=24 # reference grid remains constant across native-grid sensitivities
    target,_=target_points('medium');theta=np.linspace(0,2*np.pi,120,endpoint=False)
    tr=QDTracker(target,40000);fef.simulate_topology=tr.simulation
    rng=np.random.default_rng(seed)
    # Reuse published mutation and surrogate implementation; no model tuning.
    support=qc.Run(seed,budget=40000);support.target=target;support.theta=theta
    bank=cfg.SurrogateBank(support) if surrogate else None
    native={};plain=[]
    def key_of(desc):
        return (int(np.clip((desc[0]-qc.D1_LO)/(qc.D1_HI-qc.D1_LO)*grid,0,grid-1)),
                int(np.clip((desc[1]-qc.D2_LO)/(qc.D2_HI-qc.D2_LO)*grid,0,grid-1)))
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            if baseline:
                # Same exact total and fixed family quotas; no family is dropped.
                for i,fam in enumerate(FAMILIES):
                    tr.cap=(i+1)*5000
                    try:
                        while tr.count<tr.cap:fef.search_topology(fam,target,theta)
                    except BudgetStop:pass
            else:
                # Identical 24k search / 16k refinement split in all four arms.
                while tr.count<24000:
                    if tr.count<4000 or rng.random()<.05 or (not native and not plain):
                        fam=FAMILIES[int(rng.integers(8))];s=fef.PARAMS[fam];x=rng.uniform(s['lb'],s['ub'])
                    else:
                        parents=list(native.values()) if map_elites else plain
                        parent=parents[int(rng.integers(len(parents)))];fam=parent['topology']
                        if bank is not None and rng.random()<.5:x=bank.guided_mutation(fam,parent['x'])
                        else:x=qc.mutate(support,fam,parent['x'])
                    ev=fef.evaluate_candidate(fam,x,target,theta)
                    if ev['valid']:
                        rec={'topology':fam,'x':ev['x'].copy(),'rms_mm':ev['rms_mm']}
                        if bank is not None:bank.add(ev)
                        if map_elites:
                            key=key_of(ev['_d']);inc=native.get(key)
                            if inc is None or rec['rms_mm']<inc['rms_mm']:native[key]=rec
                        else:
                            plain.append(rec)
                            if len(plain)>576:plain=sorted(plain,key=lambda r:r['rms_mm'])[:576]
                candidates=list(native.values()) if map_elites else plain
                # Fixed family quotas in common refinement, independent of winners.
                for i,fam in enumerate(FAMILIES):
                    tr.cap=24000+(i+1)*2000
                    options=sorted([r for r in candidates if r['topology']==fam],key=lambda r:r['rms_mm'])
                    if not options:
                        s=fef.PARAMS[fam];options=[{'x':rng.uniform(s['lb'],s['ub'])}]
                    k=0
                    try:
                        while tr.count<tr.cap:
                            fef.refine_candidate(fam,options[k%len(options)]['x'],target,theta);k+=1
                    except BudgetStop:pass
    finally:fef.simulate_topology=re.simulation
    assert tr.count==40000
    strict=[]
    for r in tr.strict_archive.values():
        r['physical']=re.physical_check(r['topology'],r['x'],1920);strict.append(r)
    best_by_family=[]
    for fam in FAMILIES:
        rr=[r for r in tr.reference.values() if r['topology']==fam]
        if rr:
            r=min(rr,key=lambda r:r['rms_mm']);r['physical']=re.physical_check(fam,r['x'],1920);best_by_family.append(r)
    d={'seed':seed,'arm':label,'native_grid':grid,'reference_grid':24,'simulations':tr.count,
       'native_map':map_elites,'surrogate':surrogate,'DE':False,'shared_archive':True,
       'source_hash':sha,'strict_vectors':strict,'best_by_family':best_by_family,
       'any_strict':bool(strict),'any_strict_motion':any(r['physical']['motion_valid'] for r in strict),
       'any_quality30':any(r['physical']['quality30'] for r in strict),
       'coverage_reference':len(tr.reference)/576,
       'qd_score_reference':sum((5-min(r['rms_mm'],5))/5 for r in tr.reference.values()),
       'best_rms_mm':tr.best['rms_mm'] if tr.best else None,'seconds':time.time()-start}
    def encode(v):
        if isinstance(v,np.ndarray):return v.tolist()
        if isinstance(v,np.generic):return v.item()
        raise TypeError(type(v).__name__)
    temp=path.with_suffix('.tmp');temp.write_text(json.dumps(d,indent=2,default=encode,allow_nan=False));temp.replace(path)
    return {'name':path.name,'seconds':d['seconds'],'strict':d['any_strict'],'motion':d['any_strict_motion']}

if __name__=='__main__':print(run_qd((int(sys.argv[1]),bool(int(sys.argv[2])),bool(int(sys.argv[3])),int(sys.argv[4]),bool(int(sys.argv[5])))))
