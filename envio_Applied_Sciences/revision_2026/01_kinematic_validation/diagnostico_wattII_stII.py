# -*- coding: utf-8 -*-
"""Diagnóstico de la falla sistemática de Watt II / Stephenson II.

Hipótesis a contrastar con evidencia:
  H1 (escala): las curvas trazadoras factibles bajo los límites compactos
     tienen un tamaño característicamente distinto de 90 mm; la alineación
     rígida (sin escala) tiene un piso de error. Si el error con escala
     (similitud) colapsa a ~0, la parametrización genera la forma pero no
     el tamaño dentro de los límites.
  H2 (forma): ni siquiera con escala la familia aproxima el objetivo
     (problema de inversión/asignación del trazador).
  H3 (tamaño de muestra del pool archivado): el pool aleatorio no cubre
     la región factible de bajo error.

Salida: JSON + CSV en esta carpeta.
"""
from __future__ import annotations
import json, math, os, sys
from pathlib import Path

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, r"C:\Users\claudia.parra\Documents\Materias_romulo\investigacion\FEF-Graph\codigos")

os.environ["FEF_GEOMETRY_PROFILE"] = "compact"
os.environ["FEF_MAX_PRINTABLE_LINK_MM"] = "180.0"
os.environ["FEF_PHASE_STRIDE"] = "4"

import numpy as np
import fef_mechanism_family_synthesis as fef

OUT = AQUI
N_SAMPLE = 400
FAMILIAS = ["wattI", "wattII", "stephensonI", "stephensonII"]
TARGETS = ["easy", "medium"]


def rms_max(P, Y):
    D = P - Y
    e = np.sqrt((D ** 2).sum(axis=1))
    return float(np.sqrt((e ** 2).mean())), float(e.max())


def align_similarity(P, Y):
    """Alineación con escala (Procrustes de similitud) + fase cíclica."""
    best = (np.inf, np.inf)
    n = len(P)
    for s in range(0, n, 4):
        Ys = np.roll(Y, s, axis=0)
        Pc, Yc = P - P.mean(0), Ys - Ys.mean(0)
        H = Pc.T @ Yc
        U, sv, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        if np.linalg.det(R) < 0:
            Vt[-1] *= -1
            R = Vt.T @ U.T
        scale = float((sv.sum()) / (Pc ** 2).sum())
        Pa = scale * (R @ P.T).T + (Ys.mean(0) - scale * (R @ P.mean(0)))
        r, m = rms_max(Pa, Ys)
        if r < best[0]:
            best = (r, m)
    return best


def main():
    theta = fef.theta_samples(fef.N_POINTS)
    rng = np.random.default_rng(20261014)
    resultados = {}
    filas = []
    for tname in TARGETS:
        target, meta = fef.generate_independent_benchmark_target(theta, tname)
        for fam in FAMILIAS:
            lb, ub = fef.PARAMS[fam]["lb"], fef.PARAMS[fam]["ub"]
            anchos, rms_rig, rms_sim, valids = [], [], [], 0
            X = lb + (ub - lb) * rng.random((N_SAMPLE, len(lb)))
            for x in X:
                sim = fef.simulate_topology(fam, x, theta)
                if not sim.success:
                    continue
                valids += 1
                P = fef.curve_for_error(sim.P, len(target))
                anchos.append(fef.path_width(P))
                a = fef.best_alignment(P, target, allow_phase=True)
                rms_rig.append(a["rms"])
                rms_sim.append(align_similarity(P, target)[0])
                filas.append({"target": tname, "family": fam, "width": anchos[-1],
                              "rms_rigid": a["rms"], "rms_similarity": rms_sim[-1],
                              "x": json.dumps([float(v) for v in x])})
            anchos, rms_rig, rms_sim = map(np.asarray, (anchos, rms_rig, rms_sim))
            resultados[f"{tname}/{fam}"] = {
                "muestras_validas": int(valids),
                "ancho_curva_mediana_mm": float(np.median(anchos)),
                "ancho_curva_p05_p95_mm": [float(np.percentile(anchos, 5)), float(np.percentile(anchos, 95))],
                "rms_rigido_min_mediana": [float(rms_rig.min()), float(np.median(rms_rig))],
                "rms_similitud_min_mediana": [float(rms_sim.min()), float(np.median(rms_sim))],
                "pct_rigido_menor_1mm": float((rms_rig < 1.0).mean() * 100),
                "pct_similitud_menor_1mm": float((rms_sim < 1.0).mean() * 100),
            }
            print(f"{tname}/{fam}: validas={valids} ancho_med={np.median(anchos):.1f} mm "
                  f"rms_rigido min/med={rms_rig.min():.3f}/{np.median(rms_rig):.3f} "
                  f"rms_sim min/med={rms_sim.min():.3f}/{np.median(rms_sim):.3f}")
    (OUT / "diagnostico_wattII_stII.json").write_text(json.dumps(resultados, indent=2), encoding="utf-8")
    import csv
    with open(OUT / "muestra_aleatoria_pool.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["target", "family", "width", "rms_rigid", "rms_similarity", "x"])
        w.writeheader()
        w.writerows(filas)
    print("OK ->", OUT)


if __name__ == "__main__":
    main()
