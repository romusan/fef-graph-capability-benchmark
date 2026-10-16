"""Three worker resumable numerical campaign. Results written atomically.

Pilot is infrastructure only (seeds 871+). Formal seeds 1..30 are disjoint.
Status is saved after every completion. Failure never silently drops a seed.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1';os.environ['NUMBA_NUM_THREADS']='1'
import sys,json,time,traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE/'source'))
from campaign import run_one,FAMILIES

def jobs(stage):
    if stage=='pilot':
        return [(fam,871,'compact',target,'hybrid',2000) for target in ['medium','medium_hard'] for fam in FAMILIES]
    if stage=='factorial':
        return [(fam,s,p,t,'hybrid',40000) for t in ['medium','medium_hard']
          for p in ['compact','relaxed'] for s in range(1,31) for fam in FAMILIES]
    if stage=='new_targets':
        return [(fam,s,'compact',t,'hybrid',20000) for t in ['bernoulli','rounded_rectangle','asymmetric_dent']
          for s in range(1,31) for fam in FAMILIES]
    if stage=='easy':
        return [(fam,s,'compact','easy','hybrid',20000) for s in range(1,31) for fam in FAMILIES]
    if stage=='repeatability':
        return [('wattI',s,'compact','medium','hybrid',20000) for s in list(range(1001,1031))+list(range(2001,2031))]
    if stage=='independent':
        return [(fam,s,'compact',t,'cma',20000) for t in ['medium','medium_hard']
          for s in range(1,31) for fam in FAMILIES]
    if stage=='singularity':
        return [(fam,s,'compact','medium','hybrid_nosing',20000)
          for s in range(1,31) for fam in FAMILIES]
    raise ValueError(stage)

def main():
    stage=sys.argv[1];tasks=jobs(stage);start=time.time()
    out=HERE/('pilot_v2' if stage=='pilot' else 'repeatability_results' if stage=='repeatability' else 'results')
    status={'stage':stage,'total':len(tasks),'completed':0,'failures':[],'started_epoch':start}
    log=HERE/(stage+'.jsonl')
    with ProcessPoolExecutor(max_workers=3) as pool:
        pending={pool.submit(run_one,j,out):j for j in tasks}
        for future in as_completed(pending):
            j=pending[future]
            try:
                r=future.result();status['completed']+=1
                with log.open('a',encoding='utf8') as f:f.write(json.dumps(r)+'\n')
                print(f"{stage}: {status['completed']}/{len(tasks)} {r}",flush=True)
            except Exception:
                err={'job':j,'error':traceback.format_exc()};status['failures'].append(err)
                print(json.dumps(err),flush=True)
            status['elapsed_seconds']=time.time()-start
            temp=HERE/(stage+'_status.tmp');temp.write_text(json.dumps(status,indent=2));temp.replace(HERE/(stage+'_status.json'))
    if status['failures']:raise RuntimeError(f"{len(status['failures'])} failed tasks; repair and resume")

if __name__=='__main__':main()
