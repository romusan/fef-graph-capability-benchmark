"""Body-joint graph audit reconstructed from the implemented loop equations.

Multiple hinges at coincident coordinates are counted as independent R pairs.
Tracer offsets are points on rigid bodies, never extra links or joints.
"""
import json
from pathlib import Path

GRAPHS={
 'wattI':{'display':'CD-G','body_names':['ground O2-O4-G','crank O2-A','coupler A-B-C','rocker B-O4','dyad C-D','dyad D-G'],
    'edges':[(0,1),(1,2),(2,3),(3,0),(2,4),(4,5),(5,0)],'coincident_hinges':False},
 'wattII':{'display':'RC-G','body_names':['ground O2-O4-G','crank O2-A','coupler A-B','rocker B-O4-C','dyad C-D','dyad D-G'],
    'edges':[(0,1),(1,2),(2,3),(3,0),(3,4),(4,5),(5,0)],'coincident_hinges':False},
 'stephensonI':{'display':'A-O4','body_names':['ground O2-O4-O4','crank O2-A-A','coupler A-B','rocker B-O4','dyad A-D','dyad D-O4'],
    'edges':[(0,1),(1,2),(2,3),(3,0),(1,4),(4,5),(5,0)],'coincident_hinges':True},
 'stephensonII':{'display':'B-O2','body_names':['ground O2-O4-O2','crank O2-A','coupler A-B-B','rocker B-O4','dyad B-D','dyad D-O2'],
    'edges':[(0,1),(1,2),(2,3),(3,0),(2,4),(4,5),(5,0)],'coincident_hinges':True},
 'stephensonIII':{'display':'C-O4','body_names':['ground O2-O4-O4','crank O2-A','coupler A-B-C','rocker B-O4','dyad C-D','dyad D-O4'],
    'edges':[(0,1),(1,2),(2,3),(3,0),(2,4),(4,5),(5,0)],'coincident_hinges':True}}

def main():
    for name,d in GRAPHS.items():
        deg=[sum(i in e for e in d['edges']) for i in range(6)]
        ternary=[i for i,v in enumerate(deg) if v==3]
        assert sorted(deg)==[2,2,2,2,3,3]
        adjacent=tuple(ternary) in d['edges'] or tuple(reversed(ternary)) in d['edges']
        d.update({'degree_sequence':deg,'ternary_bodies':ternary,'ternary_adjacent':adjacent,
                  'chain_class':'Watt' if adjacent else 'Stephenson',
                  'generic_mobility':3*(6-1)-2*7})
    doc={'source':'https://www.sjf.stuba.sk/buxus/docs/pracoviska/ATC/PDF/TM_III/12_L_Synthesis.pdf',
         'source2':'https://www.sciencedirect.com/science/article/pii/S0094114X15002384',
         'graphs':GRAPHS,
         'conclusion':'Legacy IDs do not represent five distinct classical six-bar inversions. Use explicit model labels.',
         'scope':'Graph classification; coincident-joint specializations and tracer dependence stated separately.'}
    out=Path(__file__).with_name('topology_audit.json');out.write_text(json.dumps(doc,indent=2))
    print(json.dumps(doc,indent=2))

if __name__=='__main__':main()
