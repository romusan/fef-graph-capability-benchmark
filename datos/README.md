# Campaign data

Per-seed outputs of every campaign reported in the manuscript. Directory names
say which experiment each one is; the engine wrote them under timestamps, and
`MANIFIESTO.json` records the original path, the file count and a SHA-256
fingerprint of each campaign.

| Directory | Feeds | Seeds |
|---|---|---|
| `benchmark/easy_all_families` | Table 5, easy target, eight families | 30 |
| `benchmark/medium_wattI` | Tables 4 and 5, medium target, compact geometry | 30 |
| `benchmark/medium_hard_relaxed` | Table 7, medium-hard under the relaxed 260 mm envelope | 30 |
| `packaging/medium_hard_compact_wattI` | Table 7, same target under the compact 180 mm envelope | 30 |
| `packaging/medium_hard_geared_five_bar` | Table 7, geared five-bar under both envelopes | 30 per cell |
| `policy/all_families_medium` | Table 9, benchmark policy over all eight families | 30 |
| `policy/qd_arms` | Table 10, baseline and the two extension arms | 15 each |

## What each run directory contains

- `run_summary.json` — the full recorded configuration of that run: sampling
  budget, refinement depth, geometry profile, link limit, differential-evolution
  parameters, seed. This is what makes a campaign comparison checkable: two
  campaigns differ only where their summaries differ.
- `candidate_archive.csv` — the refined candidates retained per family. The
  manuscript tables select from here by lowest RMS error, the rule the published
  statistical runner uses.
- `all_valid_evaluations.csv` — every candidate that closed its loops.
- `best_accuracy.json`, `best_printable.json` — the selected designs.
- `target_path.csv`, `target_metadata.json` — the target curve. Identical across
  seeds of the same target, which is how the packaging comparison was verified
  to change the envelope and nothing else.

The QD arms in `policy/qd_arms` are flat JSON files instead: `Estrict_seed<N>`
is the baseline, `Bplus_seed<N>` the natural-budget extension and
`Bmatched_seed<N>` the approximately budget-matched one, each recording the
simulation count, the strict designs found and their parameter vectors.

## Regenerating the tables

The generators in `../tablas/` read these directories directly. No number in
the manuscript was entered by hand.

## Rebuilding this directory

```
py utilidades/empaqueta_datos.py <path-to-FEF-Graph>
```

It copies from the working tree and recomputes the manifest. The fingerprints
let you confirm that what is published here is what the tables were computed
from.
