"""B+ con PRESUPUESTO IGUALADO a la linea base.

Identico a run_bplus.py salvo el presupuesto de busqueda, reducido de
60.000 a 34.000 evaluaciones para que el total de simulaciones (busqueda
+ refinado) iguale la media medida del brazo Estrict. Responde a la
objecion de que la extension solo gana porque gasta mas.

Escribe results/Bmatched_seed<N>.json; NO toca Bplus_seed<N>.json.

Objetivo declarado ANTES de ejecutar: el paper identifica su propia brecha ---
para soldadura/mecanizado manda el criterio ESTRICTO (RMS<=1 mm Y MAX<=2 mm),
y en el objetivo medium/compact Watt I solo lo pasa en 7/30 semillas (23%),
con el resto de familias en 0/30. La meta de B+ es pasar el criterio estricto
de forma fiable Y con VARIOS disenos distintos por corrida (portafolio), no un
unico ganador.

Protocolo (identico al del paper donde aplica):
  * FEF_PHASE_STRIDE=1, LOCAL_MAXITER=100, LS_MAX_NFEV=300 (los valores del
    runner oficial para medium), objetivo medium, perfil compact.
  * Fase 1: busqueda B (MAP-Elites + mutacion guiada por surrogate), 60 000
    evaluaciones.
  * Fase 2: seleccion de hasta 16 elites diversos del archivo nativo (mejores
    2 celdas por familia presente, completado con los mejores globales) y
    refinado de cada uno con fef.refine_candidate --- el MISMO refinador del
    pipeline publicado, sin modificar.
  * Toda simulacion de las dos fases se cuenta (envoltorio en memoria del
    simulador); el total se reporta.

Metricas: mejor RMS y MAX finales, numero de disenos refinados que pasan el
criterio estricto, numero de celdas de referencia DISTINTAS entre ellos (esa
es la novedad frente a E: diversidad bajo tolerancia), y familias que pasan.

Uso:  python run_bplus.py <seed>   (semillas de evaluacion: 101, 202, 303)
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

# protocolo del paper ANTES de importar el motor
os.environ["FEF_PHASE_STRIDE"] = "1"
os.environ["FEF_LOCAL_MAXITER"] = "100"
os.environ["FEF_LS_MAX_NFEV"] = "300"
os.environ["FEF_RUN_LS_POLISH"] = "1"
# Presupuesto igualado a la media medida de la linea base (52,421
# simulaciones por semilla) menos el coste medio del refinado (18,632).
# Debe fijarse antes de importar qd_core, que lee BUDGET al importarse.
os.environ.setdefault("QD_BUDGET", "34000")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import numpy as np

import qd_core
from qd_core import BUDGET, cell_of
import qd_configs as cfg
import fef_mechanism_family_synthesis as fef

OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
N_REFINE = 16


def select_elites(arch, n=N_REFINE):
    """Mejores 2 celdas por familia presente + mejores globales hasta n."""
    per_fam: dict[str, list] = {}
    for rec in arch.cells.values():
        per_fam.setdefault(rec["topology"], []).append(rec)
    chosen, seen = [], set()
    for fam, recs in per_fam.items():
        for rec in sorted(recs, key=lambda r: r["rms_mm"])[:2]:
            chosen.append(rec)
            seen.add(id(rec))
    resto = sorted((r for r in arch.cells.values() if id(r) not in seen),
                   key=lambda r: r["rms_mm"])
    chosen.extend(resto[: max(0, n - len(chosen))])
    return sorted(chosen, key=lambda r: r["rms_mm"])[:n]


def main():
    seed = int(sys.argv[1])
    t0 = time.perf_counter()

    sims = {"n": 0}
    orig_sim = fef.simulate_topology

    def counting_sim(topo, x, theta):
        sims["n"] += 1
        return orig_sim(topo, x, theta)

    fef.simulate_topology = counting_sim
    try:
        # ---- fase 1: busqueda B ----
        run = cfg.run_map_elites(seed, guided=True)
        evals_search = run.evals
        arch = run.native_archive

        # ---- fase 2: refinado verificado ----
        elites = select_elites(arch)
        refined = []
        for k, rec in enumerate(elites, 1):
            ref = fef.refine_candidate(rec["topology"], np.asarray(rec["x"]),
                                       run.target, run.theta)
            d1, d2 = qd_core.descriptors(ref)
            refined.append({
                "topology": ref["topology"],
                "rms_mm": float(ref["rms_mm"]),
                "max_mm": float(ref["max_mm"]),
                "strict": bool(ref["rms_mm"] <= 1.0 and ref["max_mm"] <= 2.0),
                "rms1": bool(ref["rms_mm"] <= 1.0),
                "cell": cell_of(d1, d2),
                "x": np.asarray(ref["x"], float).tolist(),
                "from_rms": float(rec["rms_mm"]),
            })
            print(f"  refinado {k:2d}/{len(elites)} {rec['topology']:14s} "
                  f"{rec['rms_mm']:.3f} -> rms={ref['rms_mm']:.3f} "
                  f"max={ref['max_mm']:.3f}"
                  f"{'  ESTRICTO' if refined[-1]['strict'] else ''}",
                  flush=True)
    finally:
        fef.simulate_topology = orig_sim

    strict = [r for r in refined if r["strict"]]
    rms1 = [r for r in refined if r["rms1"]]
    best = min(refined, key=lambda r: r["rms_mm"])
    rec = {
        "config": "Bmatched", "seed": seed,
        "evals_search": evals_search,
        "sims_total": sims["n"],
        "n_refined": len(refined),
        "n_strict": len(strict),
        "n_rms1": len(rms1),
        "strict_cells": len({tuple(r["cell"]) for r in strict}),
        "strict_families": sorted({r["topology"] for r in strict}),
        "rms1_families": sorted({r["topology"] for r in rms1}),
        "best_rms_mm": best["rms_mm"], "best_max_mm": best["max_mm"],
        "best_family": best["topology"],
        "coverage_ref": run.ref.coverage(), "qd_score_ref": run.ref.qd_score(),
        "elapsed_s": time.perf_counter() - t0,
        "refined": refined,
    }
    p = OUT / f"Bmatched_seed{seed}.json"
    p.write_text(json.dumps(rec, indent=2), encoding="utf-8")
    print(json.dumps({k: rec[k] for k in
                      ("seed", "sims_total", "n_strict", "strict_cells",
                       "strict_families", "best_rms_mm", "best_max_mm")}))
    print(f"-> {p}")


if __name__ == "__main__":
    main()
