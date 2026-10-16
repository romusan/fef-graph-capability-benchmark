# FEF-Graph: corrected Applied Sciences revision

This branch supports **Planar Path Synthesis Across Eight Kinematic Models: Search Budget, Closed-Cycle Validity and Quality Diversity**, manuscript applsci-4642776.

The historical `main` branch is preserved. This revision repairs slider-crank rigid geometry and geared-five-bar one-turn periodicity, reconstructs model identities, and separates fitting accuracy, sampled cycle validity and all-dyad transmission quality. Legacy names are software identifiers, not eight distinct classical mechanism families. Historical data in `datos/` are retained as historical evidence.

## Read the paper and evidence

- [Clean revised manuscript](envio_Applied_Sciences/revision_2026/08_manuscript/manuscrito_applied_sciences.pdf)
- [Revision-highlighted manuscript](envio_Applied_Sciences/revision_2026/08_manuscript/manuscrito_applied_sciences_CAMBIOS.pdf)
- [Evidence index](envio_Applied_Sciences/revision_2026/15_auditoria_reenvio/INDEX_EVIDENCIAS.md)
- [Complete analysis](envio_Applied_Sciences/revision_2026/13_cierre_editorial/results/analysis_v2.json)
- [Statistical contrasts](envio_Applied_Sciences/revision_2026/13_cierre_editorial/results/extended_analysis_v2.json)
- [Model and repetition diagnosis](envio_Applied_Sciences/revision_2026/15_auditoria_reenvio/model_diagnostics.json)

## Experiment inventory

| Campaign | Runs | Simulator calls per run |
|---|---:|---:|
| Eight models × two targets × two profiles × 30 seeds | 960 | 40000 |
| Independent CMA-ES: eight models × two targets × 30 seeds | 480 | 20000 |
| Easy + three additional analytic paths + NACA0012 reference | 1200 | 20000 |
| Penalty ablation | 240 | 20000 |
| **Main model-level total** | **2880** | |
| Separate repetition batches | 60 | 20000 |
| Separate own-curve recovery diagnostics | 8 | 80000 |
| Separate matched-budget QD component/grid runs | 240 | 40000 total across models |

Factorial checkpoints are dependent cumulative observations, not extra independent runs. Pilot seeds do not enter formal rates. NACA0012 is published reference geometry, not a measured mechanism trajectory. All failures keep their denominators. A minimum-RMS selected vector can differ from the retained strict vector; resolution sensitivity holds the selected vector fixed.

## Verify and reproduce

Clone this branch with its directory structure intact. `python verify_publication.py` uses only the Python standard library to verify the published revision manifest, exact budgets, core/QD source fingerprints and formal counts. It does not assert acceptance of the manuscript.

The corrected numerical pipeline uses Python 3.12 and the pinned packages in `envio_Applied_Sciences/revision_2026/12_campanas_corregidas/requirements.lock.txt`. From `envio_Applied_Sciences/revision_2026/` run:

```text
python 13_cierre_editorial/summarize_v2.py
python 13_cierre_editorial/extended_analysis.py
python 13_cierre_editorial/validate_resolution_v2.py
python 13_cierre_editorial/render_tables_figures.py
python 13_cierre_editorial/finish_diagnostics.py
python 15_auditoria_reenvio/diagnose_models.py
python 15_auditoria_reenvio/envelope_diagnostics.py
python 15_auditoria_reenvio/editorial_tables.py
```

No new optimization is required to regenerate the analysis. Full run controllers and protocols are in `12_campanas_corregidas/`. The corrected Numba engine and original Python snapshot are included; the optional original C++ engine is unnecessary for this pipeline. Existing results reject incompatible source hashes. A new numerical replication should use a separate checkout and preserve supplied outputs as a reference.

The canonical LaTeX paper and all tables/figures are included in `08_manuscript/`; compile there with pdfLaTeX twice. Journal correspondence and the original review reports are provided separately through the editorial submission, so private editorial builders/backups are not needed to reproduce the numerical evidence or compile the published paper. The highlighted PDF compares against the preserved local pre-review Applied Sciences PDF; it marks inserted/replaced words and retains the clean PDF's text and pagination.

The existing repository MIT license covers the software and associated documentation. Third-party reference geometry retains its recorded source and provenance. This GitHub branch is a public access location and does not carry a DOI. Use a Git commit permalink for an immutable version. See `README_HISTORICAL.md` for the earlier campaign description.
