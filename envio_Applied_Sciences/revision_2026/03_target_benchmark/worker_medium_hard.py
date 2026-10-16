# -*- coding: utf-8 -*-
"""Worker: una corrida medium-hard de una familia / semilla / perfil.

Replica exactamente la configuracion de la campana original all-family
medium-hard (ver run_summary.json de
repo_codigo_nuevo/datos/benchmark/medium_hard_relaxed/paper_stats_run_*):
  n_random=900, shape_seeds=36, topk_refine=8, elite_neighborhoods=8,
  elite_neighbor_samples=6, DE restarts=2 x (maxiter 55, popsize 8),
  local_maxiter=220, ls_max_nfev=800, n_points=120, accept 1.0/2.0 mm,
  random_sampler=lhs, phase_stride=4, target=curva de Fourier
  "medium_hard" generada analiticamente (identica en todas las corridas,
  ancho 90 mm).

Perfiles: compact (cota de varilla <= 180 mm, corrigiendo el bound por
defecto de 210 mm) y relaxed (FEF_GEOMETRY_PROFILE=relaxed, cota 260 mm).

Reanudable: si run_summary.json ya existe en la carpeta de salida, sale
inmediatamente con codigo 0 ("ya completada").

Uso: python worker_medium_hard.py <familia> <semilla> <perfil>
  familia: fourbar | slidercrank | wattI | wattII | stephensonI |
           stephensonII | stephensonIII | gearedFiveBar
  perfil:  compact | relaxed
"""
from __future__ import annotations
import json, os, sys, time
from pathlib import Path

FAMILIAS = ["fourbar", "slidercrank", "wattI", "wattII",
            "stephensonI", "stephensonII", "stephensonIII", "gearedFiveBar"]
PERFILES = {"compact": "180.0", "relaxed": "260.0"}

if len(sys.argv) != 4 or sys.argv[1] not in FAMILIAS or sys.argv[2].isdigit() is False or sys.argv[3] not in PERFILES:
    print("Uso: python worker_medium_hard.py <familia> <semilla> <perfil compact|relaxed>")
    print("Familias:", ", ".join(FAMILIAS))
    sys.exit(2)

FAMILIA = sys.argv[1]
SEED = int(sys.argv[2])
PERFIL = sys.argv[3]

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, r"C:\Users\claudia.parra\Documents\Materias_romulo\investigacion\FEF-Graph\codigos")

# --- Config medium-hard (identica a la campana original) ---
os.environ["FEF_GEOMETRY_PROFILE"] = PERFIL
os.environ["FEF_MAX_PRINTABLE_LINK_MM"] = PERFILES[PERFIL]
os.environ["FEF_PHASE_STRIDE"] = "4"
os.environ["FEF_N_POINTS"] = "120"
os.environ["FEF_N_RANDOM"] = "900"
os.environ["FEF_SHAPE_SEEDS"] = "36"
os.environ["FEF_TOPK_REFINE"] = "8"
os.environ["FEF_RUN_DE"] = "1"
os.environ["FEF_DE_RESTARTS"] = "2"
os.environ["FEF_DE_MAXITER"] = "55"
os.environ["FEF_DE_POPSIZE"] = "8"
os.environ["FEF_LOCAL_MAXITER"] = "220"
os.environ["FEF_RUN_LS_POLISH"] = "1"
os.environ["FEF_LS_MAX_NFEV"] = "800"
os.environ["FEF_RANDOM_SAMPLER"] = "lhs"
os.environ["FEF_ELITE_NEIGHBORHOODS"] = "8"
os.environ["FEF_ELITE_NEIGHBOR_SAMPLES"] = "6"
os.environ["FEF_ACCEPT_RMS_MM"] = "1.0"
os.environ["FEF_ACCEPT_MAX_MM"] = "2.0"
os.environ["FEF_SEED"] = str(SEED)
# Sin figuras ni export CAD: campana numerica pura
os.environ["FEF_SAVE_FIGURES"] = "0"
os.environ["FEF_EXPORT_STL"] = "0"
os.environ["FEF_EXPORT_STEP"] = "0"
os.environ["FEF_EXPORT_SVG"] = "0"
os.environ["FEF_EXPORT_ZIP"] = "0"

import numpy as np
import fef_mechanism_family_synthesis as fef

SALIDA = AQUI / f"medium_hard_{PERFIL}" / f"{FAMILIA}_seed{SEED:03d}"
if (SALIDA / "run_summary.json").exists():
    print(f"{PERFIL}/{FAMILIA} seed {SEED}: ya completada, se omite")
    sys.exit(0)
SALIDA.mkdir(parents=True, exist_ok=True)

# --- Cota compacta: varilla del slider-crank de 210 a 180 mm ---
# (el perfil compacto no la ajusta por defecto; asi ambos perfiles respetan
#  la cota de fabricacion correspondiente)
if PERFIL == "compact":
    idx_l = fef.PARAMS["slidercrank"]["names"].index("l")
    if fef.PARAMS["slidercrank"]["ub"][idx_l] > 180.0:
        fef.PARAMS["slidercrank"]["ub"][idx_l] = 180.0

theta = fef.theta_samples(fef.N_POINTS)
target, meta = fef.generate_independent_benchmark_target(theta, "medium_hard")

t0 = time.time()
all_valid, refined = fef.search_topology(FAMILIA, target, theta)

import pandas as pd
filas = []
for ev in refined:
    filas.append({**{k: float(v) for k, v in zip(fef.PARAMS[FAMILIA]["names"], ev["x"])},
                  "rms_mm": ev["rms_mm"], "max_mm": ev["max_mm"],
                  "accepted": ev["accepted"], "min_delta": ev.get("min_delta", np.nan)})
pd.DataFrame(filas).to_csv(SALIDA / "candidate_archive.csv", index=False)

mejor = min(refined, key=lambda e: e["rms_mm"]) if refined else None
resumen = {
    "family": FAMILIA, "seed": SEED, "target_mode": "medium_hard",
    "geometry_profile": PERFIL, "max_printable_link_mm": float(PERFILES[PERFIL]),
    "config": {"n_random": 900, "shape_seeds": 36, "topk_refine": 8,
               "elite_neighborhoods": 8, "elite_neighbor_samples": 6,
               "de_restarts": 2, "de_maxiter": 55, "de_popsize": 8,
               "local_maxiter": 220, "ls_max_nfev": 800, "n_points": 120,
               "accept_rms_mm": 1.0, "accept_max_mm": 2.0,
               "random_sampler": "lhs", "phase_stride": 4},
    "n_valid": len(all_valid), "n_refined": len(refined),
    "best_rms_mm": mejor["rms_mm"] if mejor else None,
    "best_max_mm": mejor["max_mm"] if mejor else None,
    "strict": bool(mejor and mejor["rms_mm"] <= 1.0 and mejor["max_mm"] <= 2.0),
    "rms_pass": bool(mejor and mejor["rms_mm"] <= 1.0),
    "elapsed_seconds": time.time() - t0,
    "python_version": sys.version,
}
(SALIDA / "run_summary.json").write_text(json.dumps(resumen, indent=2), encoding="utf-8")
print(f"{PERFIL}/{FAMILIA} seed {SEED}: best RMS={resumen['best_rms_mm']:.3f} "
      f"MAX={resumen['best_max_mm']:.3f} strict={resumen['strict']} "
      f"valid={resumen['n_valid']} t={resumen['elapsed_seconds']:.0f}s")
