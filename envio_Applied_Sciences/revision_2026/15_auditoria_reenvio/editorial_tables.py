import sys,shutil
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'13_cierre_editorial'))
from render_tables_figures import latex_table
rows=[[r'Factorial: 8 models $\times$ 2 targets $\times$ 2 profiles $\times$ 30 seeds',960,'40 000'],
 [r'Independent CMA-ES: 8 $\times$ 2 $\times$ 30',480,'20 000'],
 [r'Easy + 3 new analytic + NACA0012: 8 $\times$ 5 $\times$ 30',1200,'20 000'],
 [r'Penalty ablation without singularity term: 8 $\times$ 30',240,'20 000'],
 [r'\textbf{Main model-level total}',r'\textbf{2880}',''],
 ['Additional repetition batches (separate)',60,'20 000'],
 ['Own-curve diagnostic searches (separate)',8,'80 000'],
 ['Multi-model QD policy/grid runs (separate)',240,'40 000 total']]
latex_table(ROOT/'08_manuscript/tables/v3_campaigns.tex',
 'Campaign inventory. The factorial checkpoints reuse the same searches and are not additional independent runs. The 1200-target block comprises five separate 240-run campaigns. Every numerical budget counts simulator invocations, including invalid states and refinement.',
 'tab:v3_campaigns','Campaign & Runs & Calls per run',rows,'@{}p{9.0cm}rr@{}')
shutil.copyfile(ROOT/'15_auditoria_reenvio/selected_paths_diagnosis_main.png',ROOT/'08_manuscript/figures/v3_model_diagnosis.png')
print('Editorial inventory and diagnosis figure ready')
