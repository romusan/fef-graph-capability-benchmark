"""Versioned corrected engine. Historical code/data are never modified.

Changes: rigid slider translation, integer nonzero gear ratios for a one-turn
cycle, effective rod bound, compiled equivalent geometry/alignment, counters.
"""
import os
os.environ['MPLBACKEND']='Agg'
os.environ['FEF_EXPORT_STEP']='0'
os.environ['FEF_GEOMETRY_PROFILE']='compact'
os.environ['FEF_ACCEPT_RMS_MM']='1'
os.environ['FEF_ACCEPT_MAX_MM']='2'
import copy,math
import numpy as np
import historical_engine as fef
from fast_geometry import simulate,resample,phase_alignment,aligned_metrics

FAMILIES=['fourbar','slidercrank','wattI','wattII','stephensonI','stephensonII','stephensonIII','gearedFiveBar']
BASE=copy.deepcopy(fef.PARAMS)
GEAR_RATIOS=np.array([-3.,-2.,-1.,1.,2.,3.])
SIM_CALLS=0

def configure(profile='compact'):
    fef.PARAMS=copy.deepcopy(BASE)
    fef.GEOMETRY_PROFILE=profile
    fef.MAX_PRINTABLE_LINK_MM=180. if profile=='compact' else 260.
    fef.apply_geometry_profile(fef.PARAMS)
    fef.PARAMS['slidercrank']['ub'][1]=fef.MAX_PRINTABLE_LINK_MM
    fef.PHASE_STRIDE=4
    fef.ACCEPT_RMS_MM=1.;fef.ACCEPT_MAX_MM=2.
    return fef.PARAMS

def canonical(family,x):
    spec=fef.PARAMS[family]
    x=np.clip(np.asarray(x,dtype=float),spec['lb'],spec['ub'])
    if family=='gearedFiveBar':
        x=x.copy();x[5]=GEAR_RATIOS[np.argmin(abs(GEAR_RATIOS-x[5]))]
    return x

def simulation(family,x,theta):
    global SIM_CALLS
    SIM_CALLS+=1
    ok,j,margin=simulate(FAMILIES.index(family),canonical(family,x),theta,True)
    if not ok:return fef.SimResult(False,None,None,0.,'loop closure')
    data={k:j[:,i,:] for i,k in enumerate(['A','B','C','D'])}
    if family=='slidercrank':data['O']=np.array([x[4],0.])
    else:
        data['O2']=np.zeros(2);data['O4']=np.array([x[0],0.])
    return fef.SimResult(True,j[:,4,:],data,float(margin))

def alignment(P,Y,allow_phase=True,phase_stride=4):
    r,m,s,ap=phase_alignment(P,Y,phase_stride,allow_phase)
    return {'rms':r,'max':m,'shift':s,'aligned':ap,'target_shifted':np.roll(Y,s,axis=0)}

def fixed_alignment(P,Y,shift):
    r,m,ap=aligned_metrics(P,Y,shift)
    return {'rms':r,'max':m,'shift':shift,'aligned':ap,'target_shifted':np.roll(Y,shift,axis=0)}

def evaluate_sim(family,x,sim,target,allow_phase=True):
    if not sim.success:
        return {'valid':False,'cost':1e9,'message':sim.message,'topology':family,'x':x}
    pp=resample(sim.P,len(target));a=alignment(pp,target,allow_phase)
    smooth=float(fef._smoothness_numba(a['aligned']))
    sing=max(0.,50.-sim.min_delta)/50.
    size=max(0.,float(np.sum(x[:min(6,len(x))]))-650.)/650.
    return {'valid':True,'accepted':bool(a['rms']<=1 and a['max']<=2),
       'relaxed':bool(a['rms']<=3 and a['max']<=6),
       'cost':a['rms']+.22*a['max']+.03*smooth+fef.W_SING*sing+.015*size,
       'rms_mm':a['rms'],'max_mm':a['max'],'smoothness':smooth,'min_delta':sim.min_delta,
       'shift':a['shift'],'topology':family,'x':x,'P_raw':sim.P,'P_eval':pp,
       'P_aligned':a['aligned'],'target_shifted':a['target_shifted'],'data':sim.data}

def evaluate(family,x,target,theta,allow_phase=True):
    x=canonical(family,x);sim=fef.simulate_topology(family,x,theta)
    if hasattr(sim,'revision_ev') and allow_phase:return sim.revision_ev
    return evaluate_sim(family,x,sim,target,allow_phase)

_historical_residual=fef.residual_for_topology_ls
def residual(family,x,target,theta,fixed_shift):
    ans=_historical_residual(family,x,target,theta,fixed_shift)
    if fef.W_SING==0 and ans.size==2*len(theta)+3:ans[-2]=0.
    return ans

fef.simulate_topology=simulation
fef.resample_closed_curve=resample
fef.best_alignment=alignment
fef.alignment_at_shift=fixed_alignment
fef.evaluate_candidate=evaluate
fef.residual_for_topology_ls=residual
configure()

def acute(u,v):
    den=np.linalg.norm(u,axis=-1)*np.linalg.norm(v,axis=-1)
    co=np.sum(u*v,axis=-1)/np.maximum(den,1e-300)
    return np.degrees(np.arccos(np.clip(abs(co),0,1)))

def physical_check(family,x,n=960,legacy=False):
    """Sampled full-state verification; no claim of an analytic certificate.

    Includes endpoint and every moving joint. Smooth assembly continuity is
    checked by nonzero oriented loop areas of constant sign; dense and shifted
    meshes are used by callers. Angles refer to each solved dyad's two links.
    """
    x=np.asarray(x,float);theta=np.linspace(0,2*np.pi,n+1)
    ok,j,margin=simulate(FAMILIES.index(family),x,theta,not legacy)
    if not ok:return {'n':n,'solver_valid':False,'motion_valid':False,'quality30':False,'reason':'closure_failed'}
    A,B,C,D,P=[j[:,i,:] for i in range(5)]
    endpoint=float(np.max(np.linalg.norm(j[-1]-j[0],axis=1)))
    rods=[];signed=[];angles=[]
    def rod(u,v,l):rods.append(float(np.max(abs(np.linalg.norm(u-v,axis=-1)-l))))
    def loop(u,v):
        signed.append(u[:,0]*v[:,1]-u[:,1]*v[:,0]);angles.append(float(np.min(acute(u,v))))
    if family=='slidercrank':
        O=np.array([0. if legacy else x[4],0.])
        rod(A,O,x[0]);rod(A,B,x[1])
        # Slider transmission angle is the acute angle of rod to guide normal,
        # i.e. 90 degrees minus rod obliquity to the sliding direction.
        angles.append(float(np.min(90.-acute(B-A,np.broadcast_to(np.array([1.,0.]),A.shape)))))
    elif family=='gearedFiveBar':
        O4=np.array([x[0],0.]);rod(A,np.zeros(2),x[1]);rod(B,O4,x[2]);rod(A,C,x[3]);rod(B,C,x[4]);loop(C-A,B-C)
    else:
        O4=np.array([x[0],0.]);rod(A,np.zeros(2),x[1]);rod(A,B,x[2]);rod(B,O4,x[3]);loop(B-A,O4-B)
        if family!='fourbar':
            if family in ['wattI','wattII']:G=np.array(x[6:8]);l5,l6=x[8:10]
            elif family=='stephensonIII':G=O4;l5,l6=x[6:8]
            else:G=O4 if family=='stephensonI' else np.zeros(2);l5,l6=x[4:6]
            rod(C,D,l5);rod(D,G,l6);loop(D-C,G-D)
    residual=max(rods);branch=all(np.all(s>1e-7) or np.all(s < -1e-7) for s in signed)
    motion=endpoint<=1e-6 and residual<=1e-6 and branch
    return {'n':n,'solver_valid':True,'motion_valid':bool(motion),'quality30':bool(motion and min(angles)>=30),
        'state_seam_mm':endpoint,'tracer_seam_mm':float(np.linalg.norm(P[-1]-P[0])),
        'link_residual_mm':residual,'branch_constant':bool(branch),'min_acute_angles_deg':angles,
        'min_angle_deg':min(angles),'max_joint_step_mm':float(np.max(np.linalg.norm(np.diff(j,axis=0),axis=2))),
        'min_delta':float(margin)}
