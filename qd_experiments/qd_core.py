"""Capa de calidad-diversidad (QD) sobre el motor de sintesis FEF ya verificado.

CONTEXTO IMPORTANTE (hallazgo de la Fase 1)
-------------------------------------------
El repositorio NO contiene MAP-Elites, ni archivo de nichos, ni QD-score, ni
descriptores de comportamiento. El experimento publicado es un motor de
sintesis inversa por familia (muestreo LHS + vecindarios de elites + evolucion
diferencial + refinado por minimos cuadrados) con criterio RMS<=1 mm. El
"candidate_archive.csv" es un registro plano, no un archivo MAP-Elites.

Esta capa construye ese componente QD COMO EXPERIMENTO NUEVO, reutilizando el
simulador y la funcion de evaluacion verificados sin modificarlos (solo
importa; jamas escribe en `codigos/`). Sobre esa base se prueban las tres
lineas del encargo (descriptores aprendidos, mutacion asistida por surrogate,
regla de reemplazo parametrizada) contra dos lineas base: el pipeline
existente y un MAP-Elites con descriptores definidos a mano.

DECISIONES DECLARADAS ANTES DE EJECUTAR
---------------------------------------
* Objetivo: benchmark independiente "medium", perfil geometrico compact,
  criterio del paper (RMS<=1.0 mm, MAX<=2.0 mm). Las 8 familias.
* Fitness por nicho: RMS en mm (minimizar).
* Espacio de comportamiento DE REFERENCIA (definido a mano, comun a todas las
  configuraciones para poder comparar):
    d1 = elongacion de la curva del acoplador: razon eje menor/eje mayor de la
         covarianza de P_eval (invariante a rotacion), en [0, 1].
    d2 = log10(turning_energy) de curve_complexity_metrics (complejidad de
         forma). Limites fijados con un pilotaje con semilla 999, y congelados.
* Rejilla 24 x 24. QD-score = suma sobre celdas llenas de
  (5 - min(rms, 5)) / 5, en [0, 1] por celda. Cobertura = celdas llenas / 576.
* Presupuesto: 60 000 llamadas al simulador por corrida (cada evaluacion,
  valida o no, cuesta 1). El pipeline existente corre a su presupuesto natural
  y se reporta cuantas simulaciones uso.
* Semillas de evaluacion: 101, 202, 303. Semilla de pilotaje/ajuste: 999
  (nunca se reusa para evaluar).
* Para TODAS las configuraciones, cada evaluacion se proyecta ademas en el
  archivo de referencia (archivo virtual, mejor por celda), de modo que
  cobertura y QD-score son comparables aunque la configuracion A use su propio
  espacio latente para decidir que conserva.
"""
from __future__ import annotations

import json
import os
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
CODIGOS = HERE.parent / "codigos"
sys.path.insert(0, str(CODIGOS))

# El motor verificado se configura por variables de entorno EN EL IMPORT.
os.environ.setdefault("FEF_TARGET_MODE", "medium")
os.environ.setdefault("FEF_GEOMETRY_PROFILE", "compact")
os.environ.setdefault("FEF_ACCEPT_RMS_MM", "1.0")
os.environ.setdefault("FEF_ACCEPT_MAX_MM", "2.0")
os.environ.setdefault("FEF_EXPORT_STEP", "0")
os.environ.setdefault("FEF_EXPORT_STL", "0")
os.environ.setdefault("FEF_SAVE_FIGURES", "0")
os.environ.setdefault("MPLBACKEND", "Agg")

import fef_mechanism_family_synthesis as fef  # noqa: E402

FAMILIES = list(fef.TOPOLOGIES)
GRID = 24
RMS_CAP = 5.0
BUDGET = int(os.environ.get("QD_BUDGET", "60000"))
BATCH = 64
INIT_EVALS = 4000
IMMIGRANT_P = 0.05
SIGMA = 0.12

# Limites del espacio de referencia. Fijados por el pilotaje (semilla 999,
# 6000 evaluaciones, 1547 validas) y CONGELADOS; ver pilot_bounds.json.
# d2 = log10(turning_energy).
D1_LO, D1_HI = 0.0, 1.0
D2_LO, D2_HI = -3.2, -0.6


# ------------------------------------------------------------------ soporte
def target_and_theta():
    theta = fef.theta_samples()
    target, meta = fef.generate_independent_benchmark_target(theta, "medium")
    return target, theta, meta


def descriptors(ev) -> tuple[float, float]:
    """Descriptores de referencia definidos a mano, invariantes a rotacion."""
    P = ev["P_eval"]
    C = np.cov((P - P.mean(axis=0)).T)
    w = np.linalg.eigvalsh(C)
    d1 = float(np.sqrt(max(w[0], 0.0) / max(w[1], 1e-12)))
    m = fef.curve_complexity_metrics(P)
    d2 = float(np.log10(max(m["turning_energy"], 1e-6)))
    return d1, d2


def cell_of(d1: float, d2: float) -> tuple[int, int]:
    i = int(np.clip((d1 - D1_LO) / (D1_HI - D1_LO) * GRID, 0, GRID - 1))
    j = int(np.clip((d2 - D2_LO) / (D2_HI - D2_LO) * GRID, 0, GRID - 1))
    return i, j


class ReferenceArchive:
    """Archivo virtual sobre el espacio de referencia: mejor RMS por celda.

    Toda configuracion lo alimenta con TODAS sus evaluaciones validas, de modo
    que cobertura y QD-score se miden igual para todas.
    """

    def __init__(self):
        self.cells: dict[tuple[int, int], dict] = {}

    def offer(self, ev, d1, d2):
        key = cell_of(d1, d2)
        inc = self.cells.get(key)
        if inc is None or ev["rms_mm"] < inc["rms"]:
            self.cells[key] = {"rms": ev["rms_mm"], "family": ev["topology"],
                               "x": np.asarray(ev["x"], float).copy(),
                               "d": (d1, d2)}

    def coverage(self) -> float:
        return len(self.cells) / (GRID * GRID)

    def qd_score(self) -> float:
        return float(sum((RMS_CAP - min(c["rms"], RMS_CAP)) / RMS_CAP
                         for c in self.cells.values()))

    def families_present(self) -> dict[str, int]:
        out: dict[str, int] = defaultdict(int)
        for c in self.cells.values():
            out[c["family"]] += 1
        return dict(out)


class Run:
    """Estado comun de una corrida: contador, registro y archivo de referencia."""

    def __init__(self, seed: int, budget: int = BUDGET):
        self.rng = np.random.default_rng(seed)
        fef.rng = np.random.default_rng(seed + 7)   # el motor usa su rng global
        self.target, self.theta, _ = target_and_theta()
        self.budget = budget
        self.evals = 0
        self.ref = ReferenceArchive()
        self.best = None          # mejor global por rms
        self.evals_to_pass = None  # primera eval con rms<=1 y max<=2
        self.trace = []           # (evals, coverage, qd, best_rms)

    def evaluate(self, topo: str, x: np.ndarray):
        if self.evals >= self.budget:
            return None
        self.evals += 1
        ev = fef.evaluate_candidate(topo, x, self.target, self.theta)
        if not ev["valid"]:
            return ev
        d1, d2 = descriptors(ev)
        ev["_d"] = (d1, d2)
        self.ref.offer(ev, d1, d2)
        if self.best is None or ev["rms_mm"] < self.best["rms_mm"]:
            self.best = {k: ev[k] for k in
                         ("rms_mm", "max_mm", "topology", "accepted")}
            self.best["x"] = np.asarray(ev["x"], float).tolist()
        if (self.evals_to_pass is None and ev["rms_mm"] <= 1.0
                and ev["max_mm"] <= 2.0):
            self.evals_to_pass = self.evals
        if self.evals % 2000 == 0:
            self.trace.append((self.evals, self.ref.coverage(),
                               self.ref.qd_score(),
                               self.best["rms_mm"] if self.best else None))
        return ev

    def metrics(self) -> dict:
        return {
            "evals": self.evals,
            "coverage": self.ref.coverage(),
            "qd_score": self.ref.qd_score(),
            "best_rms_mm": self.best["rms_mm"] if self.best else None,
            "best_family": self.best["topology"] if self.best else None,
            "best_accepted": bool(self.best["accepted"]) if self.best else False,
            "evals_to_pass": self.evals_to_pass,
            "cells_per_family": self.ref.families_present(),
        }


# ------------------------------------------------------- operador de variacion
def random_genome(run: Run) -> tuple[str, np.ndarray]:
    topo = FAMILIES[int(run.rng.integers(len(FAMILIES)))]
    lb = np.asarray(fef.PARAMS[topo]["lb"], float)
    ub = np.asarray(fef.PARAMS[topo]["ub"], float)
    return topo, run.rng.uniform(lb, ub)


def mutate(run: Run, topo: str, x: np.ndarray) -> np.ndarray:
    lb = np.asarray(fef.PARAMS[topo]["lb"], float)
    ub = np.asarray(fef.PARAMS[topo]["ub"], float)
    span = ub - lb
    y = x + run.rng.normal(0.0, SIGMA, size=x.shape) * span
    return np.clip(y, lb, ub)


# ------------------------------------------------- archivo MAP-Elites (nativo)
class EliteArchive:
    """Archivo MAP-Elites nativo: decide QUE se conserva y de donde se mutan
    padres. Puede vivir en el espacio de referencia (config M/B/C) o en un
    espacio latente aprendido (config A)."""

    def __init__(self, grid: int = GRID):
        self.grid = grid
        self.cells: dict[tuple[int, int], dict] = {}
        self.age: dict[tuple[int, int], int] = {}

    def keys(self):
        return list(self.cells)

    def offer(self, key, ev, replace_rule=None, neighbors_filled=0.0) -> bool:
        inc = self.cells.get(key)
        if inc is None:
            take = True
        elif replace_rule is None:
            take = ev["rms_mm"] < inc["rms_mm"]
        else:
            take = replace_rule(ev, inc, self.age.get(key, 0), neighbors_filled)
        if take:
            rec = {"rms_mm": ev["rms_mm"], "topology": ev["topology"],
                   "x": np.asarray(ev["x"], float).copy()}
            if "_sig" in ev:            # la config A necesita la firma de curva
                rec["_sig"] = ev["_sig"]
            self.cells[key] = rec
            self.age[key] = 0
        else:
            self.age[key] = self.age.get(key, 0) + 1
        return take

    def neighbors_filled_frac(self, key) -> float:
        i, j = key
        n = 0
        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                if di == dj == 0:
                    continue
                if (i + di, j + dj) in self.cells:
                    n += 1
        return n / 8.0

    def sample_parent(self, rng) -> dict:
        key = self.keys()[int(rng.integers(len(self.cells)))]
        return self.cells[key]
