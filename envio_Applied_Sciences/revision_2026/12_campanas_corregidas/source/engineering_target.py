"""Data-only target adapter for the frozen v2 search engine.

The published TMR/NASA-linked altered NACA0012 contour is an engineering
reference geometry, not a measured trajectory. No frozen engine file changes.
"""
import json,hashlib
from pathlib import Path
import numpy as np
from fast_geometry import resample
import campaign

ROOT=Path(__file__).resolve().parents[2]
RAW=ROOT/'12_campanas_corregidas/targets/n0012points_superbig_clust_fix.dat'
_original_target_points=campaign.target_points

def engineering_target(mode='naca0012_reference',n=120):
    if mode!='naca0012_reference':return _original_target_points(mode,n)
    raw=np.loadtxt(RAW,skiprows=4)[:,:2]
    assert raw.ndim==2 and len(raw)>1000 and np.isfinite(raw).all()
    assert np.linalg.norm(raw[-1]-raw[0])<1e-10
    target=resample(raw,n);target-=target.mean(axis=0);target*=90/np.ptp(target[:,0])
    return target,{'target_mode':mode,'source':'https://tmbwg.github.io/turbmodels/naca0012_val.html',
       'download_url':'https://tmbwg.github.io/turbmodels/NACA0012_grids/n0012points_superbig_clust_fix.dat',
       'source_context':'NASA-linked Turbulence Modeling Resource, altered sharp-trailing-edge NACA0012 reference geometry',
       'raw_sha256':hashlib.sha256(RAW.read_bytes()).hexdigest(),
       'target_adapter_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
       'raw_points':len(raw),'raw_closed':True,'width_mm':90.,'measured_trajectory':False,
       'preprocessing':'drop z=0, preserve contour order, closed arclength resample, center, normalize width; no smoothing or fitted mechanism'}

campaign.target_points=engineering_target

def run_engineering(job):
    family,seed=job
    out=ROOT/'12_campanas_corregidas/engineering_results';out.mkdir(exist_ok=True)
    path=out/f'naca0012_reference_compact_hybrid_{family}_seed{seed:03d}_b20000.json'
    target,meta=engineering_target()
    if path.exists():
        old=json.loads(path.read_text())['target']
        assert old['raw_sha256']==meta['raw_sha256'] and old['target_adapter_sha256']==meta['target_adapter_sha256'], 'Stale engineering target'
    return campaign.run_one((family,seed,'compact','naca0012_reference','hybrid',20000),out)

if __name__=='__main__':
    target,meta=engineering_target();print(json.dumps(meta,indent=2))
    (RAW.parent/'naca0012_provenance.json').write_text(json.dumps(meta,indent=2))
    np.savetxt(RAW.parent/'naca0012_reference_120.csv',target,delimiter=',',header='x_mm,y_mm',comments='')
