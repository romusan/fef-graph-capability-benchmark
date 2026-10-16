import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,json,time,traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE/'source'))
from qd_campaign import run_qd

def main():
    pilot='--pilot' in sys.argv
    seeds=[871] if pilot else range(1,31)
    jobs=[(s,m,g,24,False) for s in seeds for m in [False,True] for g in [False,True]]
    jobs += [(s,False,False,24,True) for s in seeds]
    if not pilot:jobs += [(s,True,True,grid,False) for s in seeds for grid in [12,18,36]]
    start=time.time();status={'total':len(jobs),'completed':0,'failures':[],'pilot':pilot}
    with ProcessPoolExecutor(max_workers=1 if pilot else 3) as pool:
        futures={pool.submit(run_qd,j):j for j in jobs}
        for fut in as_completed(futures):
            try:
                result=fut.result();status['completed']+=1;print(result,flush=True)
            except Exception:status['failures'].append({'job':futures[fut],'error':traceback.format_exc()})
            status['seconds']=time.time()-start
            p=HERE/('qd_pilot_status.json' if pilot else 'qd_status.json');p.write_text(json.dumps(status,indent=2))
    if status['failures']:raise RuntimeError(status['failures'])

if __name__=='__main__':main()
