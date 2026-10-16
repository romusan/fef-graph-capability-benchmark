# -*- coding: utf-8 -*-
"""Campaña: slider-crank con varilla <= 180 mm (perfil compacto corregido).

Motivación (obs. 1.4 del Revisor 1): la varilla del slider-crank estaba
acotada a 210 mm contra el límite compacto de 180 mm, rompiendo la
simetría entre perfiles. Esta campaña re-corre el slider-crank con la
cota de varilla fijada a 180 mm, manteniendo el resto de la
configuración idéntica a la campaña all-family medium (DE 45x9,
n_random 120, shape seeds 16, vecindades de élites 6x4, refinado top-4,
LS 300, N=120, objetivo medium generado analíticamente).

El script es reanudable: omite semillas cuya carpeta de salida ya existe.
Un proceso por semilla (misma disciplina que las campañas originales).

Uso: python corrida_slidercrank_180.py <semilla>
"""
from __future__ import annotations
import json, math, os, sys, time
from pathlib import Path

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, r"C:\Users\claudia.parra\Documents\Materias_romulo\investigacion\FEF-Graph\codigos")

# Configuración idéntica a la campaña all-family medium (run_summary.json)
os.environ["FEF_GEOMETRY_PROFILE"] = "compact"
os.environ["FEF_MAX_PRINTABLE_LINK_MM"] = "180.0"
os.environ["FEF_PHASE_STRIDE"] = "4"
os.environ["FEF_N_RANDOM"] = "120"
os.environ["FEF_SHAPE_SEEDS"] = "16"
os.environ["FEF_DE_MAXITER"] = "45"
os.environ["FEF_DE_POPSIZE"] = "9"
os.environ["FEF_TOPK_REFINE"] = "4"
os.environ["FEF_LS_MAX_NFEV"] = "300"
os.environ["FEF_LOCAL_MAXITER"] = "100"
os.environ["FEF_ELITE_NEIGHBORHOODS"] = "6"
os.environ["FEF_ELITE_NEIGHBOR_SAMPLES"] = "4"
os.environ["FEF_ACCEPT_RMS_MM"] = "1.0"
os.environ["FEF_ACCEPT_MAX_MM"] = "2.0"
os.environ["FEF_SEED"] = sys.argv[1] if len(sys.argv) > 1 else "1"

import numpy as np
import fef_mechanism_family_synthesis as fef

SEED = int(os.environ["FEF_SEED"])
SALIDA = AQUI / "slidercrank_180_medium" / f"slidercrank180_medium_seed{SEED:03d}"
if (SALIDA / "run_summary.json").exists():
    print(f"seed {SEED}: ya completada, se omite")
    sys.exit(0)
SALIDA.mkdir(parents=True, exist_ok=True)

# --- Corrección de la cota: varilla l de 210 a 180 mm ---
idx_l = fef.PARAMS["slidercrank"]["names"].index("l")
assert fef.PARAMS["slidercrank"]["ub"][idx_l] == 210.0
fef.PARAMS["slidercrank"]["ub"][idx_l] = 180.0

theta = fef.theta_samples(fef.N_POINTS)
target, meta = fef.generate_independent_benchmark_target(theta, "medium")

t0 = time.time()
all_valid, refined = fef.search_topology("slidercrank", target, theta)

import pandas as pd
filas = []
for ev in refined:
    filas.append({**{k: float(v) for k, v in zip(fef.PARAMS["slidercrank"]["names"], ev["x"])},
                  "rms_mm": ev["rms_mm"], "max_mm": ev["max_mm"],
                  "accepted": ev["accepted"], "min_delta": ev.get("min_delta", np.nan)})
df = pd.DataFrame(filas)
df.to_csv(SALIDA / "candidate_archive.csv", index=False)

mejor = min(refined, key=lambda e: e["rms_mm"]) if refined else None
resumen = {
    "seed": SEED, "topology": "slidercrank", "target_mode": "medium",
    "geometry_profile": "compact", "rod_bound_mm": 180.0,
    "config": {"de_maxiter": 45, "de_popsize": 9, "n_random": 120,
               "shape_seeds": 16, "topk_refine": 4, "ls_max_nfev": 300,
               "local_maxiter": 100, "elite_neighborhoods": 6,
               "elite_neighbor_samples": 4, "n_points": 120},
    "n_valid": len(all_valid), "n_refined": len(refined),
    "best_rms_mm": mejor["rms_mm"] if mejor else None,
    "best_max_mm": mejor["max_mm"] if mejor else None,
    "strict": bool(mejor and mejor["rms_mm"] <= 1.0 and mejor["max_mm"] <= 2.0),
    "rms_pass": bool(mejor and mejor["rms_mm"] <= 1.0),
    "elapsed_seconds": time.time() - t0,
    "python_version": sys.version,
}
(SALIDA / "run_summary.json").write_text(json.dumps(resumen, indent=2), encoding="utf-8")
print(f"seed {SEED}: best RMS={resumen['best_rms_mm']:.3f} MAX={resumen['best_max_mm']:.3f} "
      f"strict={resumen['strict']} valid={resumen['n_valid']} t={resumen['elapsed_seconds']:.0f}s")
