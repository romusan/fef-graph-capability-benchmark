"""Standard-library verification of the public scientific revision checkout."""
import json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    manifest=json.loads((ROOT/'PUBLICATION_MANIFEST_SHA256.json').read_text())
    for r in manifest['files']:
        p=ROOT/r['path'];assert p.is_file(),r['path'];assert sha(p)==r['sha256'],r['path']
    rev=ROOT/'envio_Applied_Sciences/revision_2026';camp=rev/'12_campanas_corregidas'
    names=['historical_engine.py','fast_geometry.py','revision_engine.py','campaign.py']
    core=hashlib.sha256(b''.join((camp/'source'/n).read_bytes() for n in names)).hexdigest()
    assert core==manifest['core_sha256']
    runs=list((camp/'results').glob('*.json'))+list((camp/'engineering_results').glob('*.json'))
    assert len(runs)==2880
    seen=set()
    for p in runs:
        d=json.loads(p.read_text());assert d['source_hash']==core
        assert d['simulations']==d['limit']
        key=(d['family'],d['seed'],d['profile'],d['target']['target_mode'],d['method'],d['limit'])
        assert key not in seen;seen.add(key)
    qd_source=[camp/'source'/n for n in names+['qd_campaign.py']]
    qd_source += [ROOT/'qd_experiments'/n for n in ['qd_core.py','qd_configs.py']]
    qhash=hashlib.sha256(b''.join(p.read_bytes() for p in qd_source)).hexdigest()
    qd=[json.loads(p.read_text()) for p in (camp/'qd_results').glob('*.json') if int(p.stem.split('seed')[-1]) in range(1,31)]
    assert len(qd)==240
    assert all(d['simulations']==40000 and d['source_hash']==qhash for d in qd)
    assert len(list((camp/'repeatability_results').glob('*.json')))==60
    assert len(list((rev/'11_validacion_completa/self_recovery_v2').glob('*.json')))==8
    analysis=json.loads((rev/'13_cierre_editorial/results/analysis_v2.json').read_text())
    assert analysis['complete'] and analysis['main_runs']==2880 and analysis['qd_runs']==240
    print(json.dumps({'passed':True,'verified_revision_files':len(manifest['files']),'main_runs':2880,
      'qd_runs':240,'additional_repetition_runs':60,'own_curve_cases':8,'core_sha256':core,'qd_sha256':qhash},indent=2))
if __name__=='__main__':main()
