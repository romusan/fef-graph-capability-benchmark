import json,time
from pathlib import Path
from campaign import run_one,ROOT

def main():
    records=[]
    for fam in ['fourbar','wattI','gearedFiveBar']:
        results=[]
        for rep in range(3):
            folder=ROOT/'12_campanas_corregidas'/f'determinism_{rep}'
            job=(fam,877,'compact','medium','hybrid',2000)
            run_one(job,folder)
            d=json.loads(next(folder.glob('*'+fam+'_seed877_*.json')).read_text())
            results.append({k:d[k] for k in ['simulations','valid','checkpoints','bounds','source_hash']})
        assert results[0]==results[1]==results[2],fam
        records.append({'family':fam,'seed':877,'replicates':3,'bitwise_json_numeric_identical':True})
    (ROOT/'11_validacion_completa/determinism.json').write_text(json.dumps(records,indent=2))
    print(json.dumps(records,indent=2))

if __name__=='__main__':main()
