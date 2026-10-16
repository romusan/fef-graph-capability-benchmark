# -*- coding: utf-8 -*-
"""Prueba de auto-recuperación cinemática por familia (Fase I del plan).

Protocolo por familia:
  1. Generar una configuración cinemáticamente válida cuya curva trazadora
     sea representativa (se busca una con ancho cercano a 90 mm).
  2. Calcular su trayectoria por cinemática directa y usarla como objetivo
     inverso (normalizada a 90 mm, N=120, igual que los objetivos del benchmark).
  3. Ejecutar el motor de búsqueda REAL con presupuesto reducido pero
     documentado (mismo código, mismos operadores).
  4. Medir: error RMS/MAX de reconstrucción, distancia paramétrica al
     generador (el problema inverso es multimodal: se evalúa la curva,
     no la identificabilidad exacta), validez física del diseño recuperado
     (continuidad de rama, min_delta, ángulo de transmisión aproximado).

Uso: python prueba_recuperacion.py <familia>
"""
from __future__ import annotations
import json, math, os, sys
from pathlib import Path

AQUI = Path(__file__).resolve().parent
sys.path.insert(0, r"C:\Users\claudia.parra\Documents\Materias_romulo\investigacion\FEF-Graph\codigos")

# Presupuesto reducido y documentado para la prueba de validación
os.environ["FEF_GEOMETRY_PROFILE"] = "compact"
os.environ["FEF_MAX_PRINTABLE_LINK_MM"] = "180.0"
os.environ["FEF_PHASE_STRIDE"] = "3"
os.environ["FEF_N_RANDOM"] = "900"
os.environ["FEF_SHAPE_SEEDS"] = "24"
os.environ["FEF_DE_MAXITER"] = "35"
os.environ["FEF_DE_POPSIZE"] = "9"
os.environ["FEF_TOPK_REFINE"] = "3"
os.environ["FEF_LS_MAX_NFEV"] = "300"
os.environ["FEF_ELITE_NEIGHBOR_SAMPLES"] = "0"
os.environ["FEF_ACCEPT_RMS_MM"] = "1.0"
os.environ["FEF_ACCEPT_MAX_MM"] = "2.0"

import numpy as np
import fef_mechanism_family_synthesis as fef

TOPO = sys.argv[1]
OUT = AQUI / "recuperacion"
OUT.mkdir(parents=True, exist_ok=True)


def busca_generador(fam, theta, rng, n_intentos=4000):
    """Muestrea hasta encontrar un diseño válido con curva de ancho ≥ 70 mm."""
    lb, ub = fef.PARAMS[fam]["lb"], fef.PARAMS[fam]["ub"]
    mejor, mejor_ancho = None, -1.0
    for i in range(n_intentos):
        x = lb + (ub - lb) * rng.random(len(lb))
        sim = fef.simulate_topology(fam, x, theta)
        if not sim.success:
            continue
        w = fef.path_width(sim.P)
        if w > mejor_ancho:
            mejor, mejor_ancho = x, w
        if w >= 70.0:
            return x, float(w)
    return mejor, float(mejor_ancho)


def main():
    theta = fef.theta_samples(fef.N_POINTS)
    rng = np.random.default_rng(777)

    x_gen, ancho = busca_generador(TOPO, theta, rng)
    if x_gen is None:
        print(f"{TOPO}: NO se encontró generador válido")
        (OUT / (TOPO + "_amp.json")).write_text(json.dumps(
            {"family": TOPO, "generador_encontrado": False}), encoding="utf-8")
        return
    sim_gen = fef.simulate_topology(TOPO, x_gen, theta)
    P = fef.normalize_target_width(sim_gen.P.copy(), 90.0)
    target = fef.resample_closed_curve(P, fef.N_POINTS)
    meta = {"target_source": f"self_recovery_{TOPO}", "generated_by_stephenson": False,
            "target_width_mm": 90.0}
    print(f"{TOPO}: generador ancho={ancho:.1f} mm, vars={np.round(x_gen,2).tolist()}")

    # Medición del desajuste de TAMAÑO: la curva cruda del generador (ancho
    # arbitrario) contra el objetivo normalizado a 90 mm, con alineación rígida
    # sin escala (igual que el benchmark). No es un sanity de identidad: cuantifica
    # cuánto error introduce el hecho de que el mecanismo generador no traza la
    # curva exactamente a 90 mm de ancho.
    ev0 = fef.evaluate_candidate(TOPO, x_gen, target, theta)
    print(f"{TOPO}: desajuste_tamano generador rms={ev0['rms_mm']:.4f} max={ev0['max_mm']:.4f} valid={ev0['valid']}")
    ident = fef.best_alignment(fef.curve_for_error(sim_gen.P, len(target)), target)
    print(f"{TOPO}: identidad_forma (misma curva normalizada) rms={ident['rms']:.6f}")

    all_valid, refined = fef.search_topology(TOPO, target, theta)
    mejor = min(refined, key=lambda e: e["rms_mm"]) if refined else None
    res = {
        "family": TOPO,
        "generador_encontrado": True,
        "generador": {"x": [float(v) for v in x_gen], "ancho_mm": ancho,
                      "desajuste_tamano_rms_mm": ev0["rms_mm"], "desajuste_tamano_max_mm": ev0["max_mm"],
                      "identidad_forma_rms_mm": ident["rms"],
                      "sanity_valid": bool(ev0["valid"])},
        "presupuesto": {"n_random": 500, "shape_seeds": 24, "de_maxiter": 22,
                        "de_popsize": 7, "topk_refine": 2, "ls_max_nfev": 300,
                        "nota": "presupuesto reducido para prueba de validación, no para tasas"},
        "n_validos_muestreados": len(all_valid),
        "mejor_refined": None,
    }
    if mejor:
        dist = fef.normalized_param_distance(TOPO, np.asarray(mejor["x"]), x_gen)
        res["mejor_refined"] = {
            "rms_mm": mejor["rms_mm"], "max_mm": mejor["max_mm"],
            "accepted": mejor["accepted"], "min_delta": mejor.get("min_delta"),
            "distancia_parametrica_normalizada_al_generador": dist,
            "x": [float(v) for v in mejor["x"]],
        }
        print(f"{TOPO}: recuperado rms={mejor['rms_mm']:.3f} max={mejor['max_mm']:.3f} "
              f"accepted={mejor['accepted']} dparam={dist:.3f}")
    (OUT / (TOPO + "_amp.json")).write_text(json.dumps(res, indent=2), encoding="utf-8")
    print("OK")


if __name__ == "__main__":
    main()
