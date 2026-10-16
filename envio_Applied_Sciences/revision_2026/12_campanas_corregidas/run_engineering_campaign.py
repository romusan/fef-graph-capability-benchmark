import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
import sys,json,time,traceback
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
HERE=Path(__file__).resolve().parent;sys.path.insert(0,str(HERE/'source'))
from engineering_target import run_engineering
from campaign import FAMILIES

def main():
    start=time.time();status={'total':240,'completed':0,'failures':[]}
    jobs=[(f,s) for s in range(1,31) for f in FAMILIES]
    with ProcessPoolExecutor(max_workers=3) as pool:
        futures={pool.submit(run_engineering,j):j for j in jobs}
        for fut in as_completed(futures):
            try:
                r=fut.result();status['completed']+=1;print(r,flush=True)
            except Exception:status['failures'].append({'job':futures[fut],'error':traceback.format_exc()})
            status['seconds']=time.time()-start
            temp=HERE/'engineering_status.tmp';temp.write_text(json.dumps(status,indent=2));temp.replace(HERE/'engineering_status.json')
    if status['failures']:raise RuntimeError(status['failures'])

if __name__=='__main__':main()
