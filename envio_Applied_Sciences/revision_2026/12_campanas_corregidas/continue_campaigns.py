"""Dependency-driven continuation of this authorized revision, not a scheduler.

Waits for the already running factorial to finish, then runs the remaining
stages serially, each with at most three numerical workers. No publication.
"""
import json,time,subprocess,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent

def main():
    while True:
        p=HERE/'factorial_status.json'
        if p.exists():
            d=json.loads(p.read_text())
            if d['completed']+len(d['failures'])==d['total']:
                if d['failures']:raise RuntimeError('Factorial errors require repair before continuation')
                break
        time.sleep(20)
    commands=[['run_campaigns.py','easy'],['run_campaigns.py','repeatability'],['run_campaigns.py','independent'],
              ['run_campaigns.py','new_targets'],['run_engineering_campaign.py'],['run_campaigns.py','singularity'],
              ['source/self_recovery_v2.py'],['run_qd_campaigns.py']]
    for i,args in enumerate(commands):
        if args[0]=='run_qd_campaigns.py':
            pilot=json.loads((HERE/'qd_pilot_status.json').read_text())
            if pilot['completed']!=pilot['total'] or pilot['failures']:
                raise RuntimeError('QD infrastructure pilot must pass before formal seeds')
        log=HERE/f'continuation_{i}_{Path(args[0]).stem}.log'
        with log.open('w',encoding='utf8') as f:
            r=subprocess.run([sys.executable,'-X','utf8',str(HERE/args[0]),*args[1:]],stdout=f,stderr=subprocess.STDOUT,cwd=HERE)
        status={'last_command':args,'returncode':r.returncode,'completed_stages':i+(r.returncode==0),'total_stages':len(commands)}
        (HERE/'continuation_status.json').write_text(json.dumps(status,indent=2))
        if r.returncode:raise RuntimeError(f'Stage failed; inspect {log}')

if __name__=='__main__':main()
