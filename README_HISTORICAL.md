# FEF-Graph — code for the capability-benchmark experiments

Code accompanying *Which Planar Mechanism Families Can Draw a Given Closed
Curve? A Capability Benchmark under Explicit Packaging Limits*.

This repository holds **only the campaign, aggregation and verification code
written for the experiments added to the manuscript**. The synthesis engine
itself (`fef_mechanism_family_synthesis.py`, the C++ geometry core and the
MAP-Elites harness) is not duplicated here; these scripts import it and expect
to run beside it.

## What each campaign does

| Script | Experiment |
|---|---|
| `campanas/c2_qd_15_semillas.ps1` | Extends both arms of the policy comparison from 3 to 15 seeds. Runs `run_bplus.py` (QD) and `run_e_strict.py` (baseline) on seeds 404–1515. |
| `campanas/c3_presupuesto_igualado.ps1` + `run_bplus_matched.py` | Repeats the QD arm with the archive-search budget cut from 60,000 to 34,000 evaluations, so that total simulations land within 7% of the baseline. Writes `Bmatched_seed<N>.json`; never touches `Bplus_*`. |
| `campanas/c4_medium_hard_compacto.ps1` | Medium-hard target under the **compact** 180 mm envelope, Watt I, 30 seeds. Replicates the published relaxed campaign field for field, changing only the packaging profile. |
| `campanas/c5_medium_hard_por_familia.ps1` | Generalises C4 over family and profile. Used for the geared five-bar under both envelopes, completing a two-family by two-envelope design with 30 seeds per cell. |

The `*_manifiesto.json` files record the seed list and a SHA-256 fingerprint of
the source files each campaign ran against.

## Table generators

Every number in the manuscript tables comes from these, not from hand editing:

- `tablas/aggregate_all_families.py` — all eight families on the medium target
  (Table 9). Per-family candidates are selected by lowest RMS, the same rule
  `extract_family_rows()` uses in the published statistical runner.
- `tablas/make_qd_extension_table.py` — the three arms of the policy comparison
  (Table 10), with Wilson intervals, Fisher exact tests and Holm adjustment
  across the four comparisons of that table.
- `tablas/make_packaging_table.py` — the two-by-two packaging design (Table 7).
  The envelopes share seeds 1–30, so the within-family contrasts are **paired**:
  exact McNemar for the rates, Wilcoxon signed-rank for the error. Between
  families the samples are independent and Fisher is used.

## Data

`datos/` holds the per-seed outputs of every campaign reported in the paper,
19.5 MB across 1,680 files, with a SHA-256 fingerprint per campaign in
`datos/MANIFIESTO.json`. See `datos/README.md` for the layout and for which
table each campaign feeds.

## Verification

`verificacion/verifica_wattII.py` re-simulates, in an independent process, the
single Watt II design that meets the strict criterion in Table 9. It reproduces
the archived errors exactly (RMS 0.6488 mm, maximum 1.5278 mm), which is what
supports counting four strict-capable families rather than three.

Note the parameter order: Watt II takes twelve variables and the archive columns
are named `r_c` and `phi_c`. Building the vector with ten variables evaluates a
different mechanism.

## Running the campaigns on Windows

Invoke the PowerShell launchers with `-Command`, not `-File`:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -Command "& '.\campanas\c5_medium_hard_por_familia.ps1' -Topologia gearedFiveBar -Perfil compact -Semillas (2..30) -Max 3"
```

With `-File`, arguments arrive as separate strings and an array parameter does
not bind: `-Semillas 2,3` is read as the single seed 23.

Long campaigns should run detached, with `utilidades/keep_awake.ps1` alongside
them. It uses `SetThreadExecutionState`, which reverts on exit and does not
alter the machine's power plan. Without it a suspended machine silently adds
hours of wall-clock time; the runs survive suspension, but the wait does not.

## Reproducibility notes

- Campaign outputs are written to new directories per seed, so re-running never
  overwrites an archived campaign.
- All simulations, including the refinement stage, are counted and reported.
- `utilidades/aplana_para_editorial_manager.py` flattens the manuscript sources
  into a single directory, which Elsevier's submission system requires.

## Licence

MIT. See `LICENSE`. The code may be reused, modified and redistributed,
including commercially, provided the copyright notice is retained.

If you use it in academic work, please cite the accompanying paper.

