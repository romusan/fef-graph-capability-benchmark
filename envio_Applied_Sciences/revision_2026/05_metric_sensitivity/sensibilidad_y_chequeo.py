# -*- coding: utf-8 -*-
"""Sensibilidad de tolerancias (Fase V) + validación cinemática post-hoc.

A) Tolerancias: se recalculan tasas de éxito de las campañas archivadas bajo
   rejillas de tolerancia alternativas (RMS 0.5/1.0/1.5/2.0 mm;
   Emax 1/2/3/4 mm) usando los pares (rms, max) por semilla de los mejores
   diseños refinados. Datos: reposición directa de CSV archivados.

B) Validación cinemática post-hoc de los diseños reportados: se re-simula
   cada mejor diseño por semilla y se verifican
     - cierre cinemático en todo el ciclo (ya es condición de validez),
     - continuidad de rama (sin salto entre configuraciones),
     - margen de singularidad min_delta,
     - ángulo de transmisión mínimo del primer bucle tipo cuatro barras
       (solo familias cuyo primer bucle es O2-A-B-O4; no se aplica a
       slider-crank ni geared five-bar, cuya arquitectura difiere).
   Salida: JSON con tasas de paso por campaña/familia.
"""
from __future__ import annotations
import json, math, os, re, sys
from pathlib import Path
import numpy as np
import pandas as pd

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, r"C:\Users\claudia.parra\Documents\Materias_romulo\investigacion\FEF-Graph\codigos")
os.environ["FEF_GEOMETRY_PROFILE"] = "compact"
os.environ["FEF_MAX_PRINTABLE_LINK_MM"] = "180.0"
os.environ["FEF_PHASE_STRIDE"] = "4"

RAIZ = Path(r"C:\Users\claudia.parra\Documents\Materias_romulo\investigacion\FEF-Graph\repo_codigo_nuevo\datos")
SAL = AQUI

RMS_REJILLA = [0.5, 1.0, 1.5, 2.0]
MAX_REJILLA = [1.0, 2.0, 3.0, 4.0]

CAMPANAS = [
    ("easy_all8", "benchmark/easy_all_families", "paper_stats_run_easy_all_families_seed*", None),
    ("medium_inicial_wattI", "benchmark/medium_wattI", "paper_stats_run_wattI_medium_seed*", None),
    ("medium_all8", "policy/all_families_medium", "p1_medium_all8_seed*", None),
]


def mejor_por_familia(d):
    out = {}
    for fam, g in d.groupby("topology"):
        out[fam] = g.loc[g.rms_mm.idxmin()]
    return out


# ---------------- A) tolerancias ----------------
tolerancias = {}
for nombre, directorio, patron, _ in CAMPANAS:
    for d in sorted(RAIZ.glob(f"{directorio}/{patron}")):
        f = d / "candidate_archive.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f)
        for fam, row in mejor_por_familia(df).items():
            clave = f"{nombre}/{fam}"
            t = tolerancias.setdefault(clave, {"n": 0, "rms": [], "max": []})
            t["n"] += 1
            t["rms"].append(float(row.rms_mm))
            t["max"].append(float(row.max_mm))
            t.setdefault("min_delta", []).append(float(row.min_delta))

tabla_tol = {}
for clave, t in tolerancias.items():
    rms, mx = np.array(t["rms"]), np.array(t["max"])
    fila = {"n": t["n"]}
    for r in RMS_REJILLA:
        for m in MAX_REJILLA:
            fila[f"rms<={r}&max<={m}"] = int(((rms <= r) & (mx <= m)).sum())
    fila["rms_mediana"] = float(np.median(rms))
    tabla_tol[clave] = fila

with open(SAL / "sensibilidad_tolerancias.json", "w", encoding="utf-8") as f:
    json.dump(tabla_tol, f, indent=2)
print("A) tolerancias OK:", len(tabla_tol), "campaña/familia")

# ---------------- B) validación cinemática post-hoc ----------------
import fef_mechanism_family_synthesis as fef

theta = fef.theta_samples(fef.N_POINTS)
CUATRO_BAR_LOOP1 = {"fourbar", "wattI", "wattII", "stephensonI", "stephensonII", "stephensonIII"}


def transmision_min_loop1(sim_data, A, B, O2, O4):
    """Ángulo de transmisión mínimo (deg) del bucle O2-A-B-O4:
    ángulo entre la biela AB y el balancín O4B."""
    mus = []
    for a, b in zip(A, B):
        u_ab = a - b          # dirección de empuje aproximada (biela)
        u_rock = O4 - b       # balancín
        na, nr = np.linalg.norm(u_ab), np.linalg.norm(u_rock)
        if na < 1e-9 or nr < 1e-9:
            continue
        c = np.dot(u_ab, u_rock) / (na * nr)
        mus.append(math.degrees(math.acos(max(-1.0, min(1.0, c)))))
    return float(min(mus)) if mus else float("nan")


def verifica(fam, x, theta):
    sim = fef.simulate_topology(fam, x, theta)
    if not sim.success:
        return {"valido": False}
    P = sim.P
    # continuidad de rama: máximo salto entre posiciones consecutivas del trazador
    saltos = np.linalg.norm(np.diff(np.vstack([P, P[0]]), axis=0), axis=1)
    resultado = {
        "valido": True,
        "min_delta": float(sim.min_delta),
        "salto_max_mm": float(saltos.max()),
        "ancho_mm": float(fef.path_width(P)),
    }
    if fam in CUATRO_BAR_LOOP1:
        A, B = sim.data["A"], sim.data["B"]
        resultado["transmision_min_deg"] = transmision_min_loop1(None, A, B, sim.data["O2"], sim.data["O4"])
    return resultado


chequeo = {}
NAMES = {"fourbar": ["L1", "L2", "L3", "L4", "rP", "phiP"], }
for nombre, directorio, patron, _ in CAMPANAS:
    for d in sorted(RAIZ.glob(f"{directorio}/{patron}")):
        f = d / "candidate_archive.csv"
        if not f.exists():
            continue
        df = pd.read_csv(f)
        for fam, row in mejor_por_familia(df).items():
            nombres = fef.PARAMS[fam]["names"]
            if not all(c in df.columns for c in nombres):
                continue
            x = np.array([float(row[c]) for c in nombres])
            v = verifica(fam, x, theta)
            clave = f"{nombre}/{fam}"
            acum = chequeo.setdefault(clave, {"n": 0, "validos": 0, "min_delta": [],
                                              "salto_max": [], "transmision_min": []})
            acum["n"] += 1
            if v.get("valido"):
                acum["validos"] += 1
                acum["min_delta"].append(v["min_delta"])
                acum["salto_max"].append(v["salto_max_mm"])
                if "transmision_min_deg" in v:
                    acum["transmision_min"].append(v["transmision_min_deg"])

resumen_chequeo = {}
for clave, a in chequeo.items():
    r = {"n": a["n"], "validos": a["validos"],
         "min_delta_mediana": float(np.median(a["min_delta"])) if a["min_delta"] else None,
         "salto_max_mediana_mm": float(np.median(a["salto_max"])) if a["salto_max"] else None}
    if a["transmision_min"]:
        tm = np.array(a["transmision_min"])
        r["transmision_min_mediana_deg"] = float(np.median(tm))
        r["transmision_min_p5_deg"] = float(np.percentile(tm, 5))
        r["frac_transmision_ge_40deg"] = float((tm >= 40).mean())
        r["frac_transmision_ge_30deg"] = float((tm >= 30).mean())
    resumen_chequeo[clave] = r

with open(SAL / "validacion_cinematica_posthoc.json", "w", encoding="utf-8") as f:
    json.dump(resumen_chequeo, f, indent=2)
print("B) chequeo OK:", len(resumen_chequeo), "campaña/familia")
for k in ["medium_all8/wattI", "medium_all8/wattII", "medium_all8/stephensonII",
          "easy_all8/wattII", "medium_inicial_wattI/wattI"]:
    if k in resumen_chequeo:
        print(" ", k, json.dumps(resumen_chequeo[k]))
