"""Read actual lineno anchors from the compiled clean manuscript."""
import sys,re,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'13_cierre_editorial'))
from editorial_support import KEYS,MAP
source=ROOT/'08_manuscript/manuscrito_applied_sciences.tex'
aux=source.with_suffix('.aux').read_text(encoding='utf8')
found={}
for key,end,line,page in re.findall(r'\\newlabel\{rev:([^:]+):(start|end)\}\{\{(\d+)\}\{(\d+)\}',aux):
    found[(key,end)]=(int(line),int(page))
locations={}
for key in list(KEYS.values())+['ack']:
    start=found[(key,'start')];end=found[(key,'end')]
    assert start[0]<=end[0] and start[1]<=end[1],(key,start,end)
    locations[key]={'first_line':start[0],'last_line':end[0],'first_page':start[1],'last_page':end[1]}
out={'source_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),
     'pdf_sha256':hashlib.sha256(source.with_suffix('.pdf').read_bytes()).hexdigest(),
     'locations':locations,'comment_keys':MAP,'basis':'lineno labels compiled in the revised clean manuscript, not source-file line indices'}
(ROOT/'13_cierre_editorial/results/revised_line_locations.json').write_text(json.dumps(out,indent=2))
print(json.dumps(locations,indent=2))
