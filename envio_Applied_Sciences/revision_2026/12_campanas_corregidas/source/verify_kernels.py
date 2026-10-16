"""Non-tautological equivalence and mechanical regression tests."""
import json,os,time
from pathlib import Path
os.environ['FEF_EXPORT_STEP']='0';os.environ['FEF_GEOMETRY_PROFILE']='compact'
import numpy as np
import historical_engine as old
from fast_geometry import simulate,resample,phase_alignment

def main():
    rng=np.random.default_rng(871);theta=np.linspace(0,2*np.pi,120,endpoint=False)
    counts={};max_error=0.;metric_error=0.;saved=[]
    for fi,family in enumerate(old.PARAMS):
        spec=old.PARAMS[family];counts[family]=0
        for rep in range(100):
            x=rng.uniform(spec['lb'],spec['ub']);s=old.simulate_topology(family,x,theta)
            ok,j,m=simulate(fi,x,theta,False)
            assert ok==s.success,(family,rep,'validity')
            if ok:
                err=float(np.max(abs(j[:,4]-s.P)));max_error=max(max_error,err)
                assert err<1e-8,(family,err)
                assert abs(m-s.min_delta)<1e-6
                p=resample(s.P,120);po=old.resample_closed_curve(s.P,120)
                assert np.max(abs(p-po))<1e-8
                y=old.generate_independent_benchmark_target(theta,'medium')[0]
                a=old.best_alignment(p,y,phase_stride=4);r,m,shift,pa=phase_alignment(p,y,4)
                metric_error=max(metric_error,abs(r-a['rms']),abs(m-a['max']))
                assert abs(r-a['rms'])<1e-8 and abs(m-a['max'])<1e-8
                counts[family]+=1
                if family not in ['slidercrank','gearedFiveBar'] and counts[family]<=3:
                    ev=old.evaluate_candidate(family,x,y,theta)
                    saved.append((family,x.copy(),y.copy(),ev['rms_mm'],ev['max_mm'],ev['cost']))
    from revision_engine import physical_check,canonical,configure
    import revision_engine as re
    for fam,x,y,r,m,cost in saved:
        ev=re.evaluate(fam,x,y,theta)
        assert abs(ev['rms_mm']-r)<1e-8 and abs(ev['max_mm']-m)<1e-8
        assert abs(ev['cost']-cost)<1e-8
    x=np.array([30.,100.,.5,20.,15.])
    assert physical_check('slidercrank',x,120,True)['link_residual_mm']>1
    assert physical_check('slidercrank',x,120)['link_residual_mm']<1e-10
    x=np.array([100.,20.,20.,90.,90.,.93,.5,10.,0.])
    assert physical_check('gearedFiveBar',x,120,True)['state_seam_mm']>1
    assert not physical_check('gearedFiveBar',x,120,True)['motion_valid']
    assert canonical('gearedFiveBar',x)[5]==1
    configure('relaxed')
    assert re.fef.PARAMS['slidercrank']['ub'][1]==260
    rec={'valid_examples':counts,'max_point_error_mm':max_error,'max_metric_error_mm':metric_error,
         'tests':'800 random cases, legacy geometry/metrics equivalence and physical repairs',
         'passed':True}
    (Path(__file__).parents[2]/'11_validacion_completa/kernel_verification.json').write_text(json.dumps(rec,indent=2))
    print(json.dumps(rec,indent=2))

if __name__=='__main__':main()
