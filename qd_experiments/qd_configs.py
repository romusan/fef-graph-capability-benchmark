"""Configuraciones del experimento QD: E, M, A, B y C.

E  pipeline existente (verificado) instrumentado, presupuesto natural.
M  MAP-Elites con los descriptores de referencia definidos a mano.
A  MAP-Elites con descriptores APRENDIDOS: autoencoder pequeno (numpy puro,
   sin torch) sobre curvas del acoplador normalizadas; estilo AURORA con un
   re-entrenamiento y re-codificacion del archivo a mitad de presupuesto.
B  = M + mutacion asistida: mitad de las mutaciones guiadas por el gradiente
   de un surrogate local (regresion ridge ponderada sobre vecinos evaluados).
   El genotipo por familia es un vector real acotado, asi que el gradiente del
   surrogate se aplica directamente; no hace falta pasar por el latente de A.
C  = M + regla de reemplazo parametrizada (logistica sobre 4 rasgos), con los
   parametros elegidos por busqueda aleatoria en la semilla de ajuste 999 y
   CONGELADOS antes de tocar las semillas de evaluacion.
"""
from __future__ import annotations

import json
import sys
from collections import deque
from pathlib import Path

import numpy as np

from qd_core import (BATCH, BUDGET, FAMILIES, GRID, IMMIGRANT_P, INIT_EVALS,
                     EliteArchive, Run, cell_of, mutate, random_genome)
import qd_core
import fef_mechanism_family_synthesis as fef


# =====================================================================
# M — MAP-Elites, descriptores manuales
# =====================================================================
def run_map_elites(seed: int, budget: int = BUDGET, replace_rule=None,
                   guided: bool = False) -> Run:
    run = Run(seed, budget)
    arch = EliteArchive()
    sur = SurrogateBank(run) if guided else None

    while run.evals < min(INIT_EVALS, budget):
        topo, x = random_genome(run)
        ev = run.evaluate(topo, x)
        if ev and ev.get("valid"):
            _offer_native(arch, ev, replace_rule)
            if sur:
                sur.add(ev)

    while run.evals < budget:
        for _ in range(BATCH):
            if run.evals >= budget:
                break
            if not arch.cells or run.rng.random() < IMMIGRANT_P:
                topo, x = random_genome(run)
            else:
                parent = arch.sample_parent(run.rng)
                topo = parent["topology"]
                if sur and run.rng.random() < 0.5:
                    x = sur.guided_mutation(topo, parent["x"])
                else:
                    x = mutate(run, topo, parent["x"])
            ev = run.evaluate(topo, x)
            if ev and ev.get("valid"):
                _offer_native(arch, ev, replace_rule)
                if sur:
                    sur.add(ev)
    run.native_archive = arch
    return run


def _offer_native(arch: EliteArchive, ev, replace_rule):
    key = cell_of(*ev["_d"])
    arch.offer(key, ev, replace_rule=replace_rule,
               neighbors_filled=arch.neighbors_filled_frac(key))


# =====================================================================
# B — surrogate local por familia
# =====================================================================
class SurrogateBank:
    """Un surrogate local por familia: ridge ponderada sobre los k vecinos mas
    cercanos del padre en parametros normalizados. Su gradiente da la
    direccion de descenso estimada del RMS."""

    K = 30
    MIN_DATA = 120
    BUF = 4000

    def __init__(self, run: Run):
        self.run = run
        self.data = {t: deque(maxlen=self.BUF) for t in FAMILIES}
        self.bounds = {t: (np.asarray(fef.PARAMS[t]["lb"], float),
                           np.asarray(fef.PARAMS[t]["ub"], float))
                       for t in FAMILIES}

    def add(self, ev):
        lb, ub = self.bounds[ev["topology"]]
        z = (np.asarray(ev["x"], float) - lb) / (ub - lb)
        self.data[ev["topology"]].append((z, min(ev["rms_mm"], 25.0)))

    def guided_mutation(self, topo: str, x: np.ndarray) -> np.ndarray:
        lb, ub = self.bounds[topo]
        d = self.data[topo]
        if len(d) < self.MIN_DATA:
            return mutate(self.run, topo, x)
        z0 = (np.asarray(x, float) - lb) / (ub - lb)
        Z = np.asarray([p[0] for p in d])
        y = np.asarray([p[1] for p in d])
        dist = np.linalg.norm(Z - z0, axis=1)
        idx = np.argsort(dist)[: self.K]
        Zk, yk, dk = Z[idx], y[idx], dist[idx]
        w = np.exp(-(dk / (dk.mean() + 1e-9)) ** 2)
        A = np.hstack([Zk - z0, np.ones((len(idx), 1))])
        W = np.diag(w)
        # ridge ponderada: (A'WA + lam I) beta = A'W y
        lam = 1e-3
        M = A.T @ W @ A + lam * np.eye(A.shape[1])
        beta = np.linalg.solve(M, A.T @ W @ yk)
        g = beta[:-1]
        gn = np.linalg.norm(g)
        if gn < 1e-9:
            return mutate(self.run, topo, x)
        eta = float(np.exp(self.run.rng.uniform(np.log(0.03), np.log(0.3))))
        z1 = z0 - eta * g / gn + self.run.rng.normal(0, 0.02, size=z0.shape)
        return np.clip(z1, 0, 1) * (ub - lb) + lb


# =====================================================================
# A — descriptores aprendidos (autoencoder numpy, estilo AURORA)
# =====================================================================
N_CURVE = 32   # puntos por curva; entrada del AE = 64 valores


def curve_signature(P: np.ndarray, rng=None) -> np.ndarray:
    """Curva del acoplador normalizada: remuestreada, centrada, escala unidad
    y orientada por PCA (invariante a rotacion y a escala)."""
    Q = fef.resample_closed_curve(P, N_CURVE)
    Q = Q - Q.mean(axis=0)
    s = np.sqrt((Q ** 2).sum(axis=1).mean()) + 1e-12
    Q = Q / s
    C = np.cov(Q.T)
    _, V = np.linalg.eigh(C)
    R = V[:, ::-1]
    if np.linalg.det(R) < 0:
        R[:, 1] = -R[:, 1]
    Q = Q @ R
    # fase canonica: empezar en el punto de mayor x
    k = int(np.argmax(Q[:, 0]))
    Q = np.roll(Q, -k, axis=0)
    return Q.reshape(-1)


class TinyAE:
    """Autoencoder 64-32-2-32-64 (tanh) en numpy puro con Adam.

    Sin torch/sklearn en el entorno verificado, un AE minimo propio es la
    opcion honesta; PCA seria el caso lineal del mismo enfoque.
    """

    def __init__(self, d_in=2 * N_CURVE, h=32, z=2, seed=0):
        r = np.random.default_rng(seed)
        s = 1.0 / np.sqrt(d_in)
        self.W1 = r.normal(0, s, (d_in, h)); self.b1 = np.zeros(h)
        self.W2 = r.normal(0, 0.2, (h, z));  self.b2 = np.zeros(z)
        self.W3 = r.normal(0, 0.5, (z, h));  self.b3 = np.zeros(h)
        self.W4 = r.normal(0, 0.2, (h, d_in)); self.b4 = np.zeros(d_in)
        self.params = [self.W1, self.b1, self.W2, self.b2,
                       self.W3, self.b3, self.W4, self.b4]
        self.m = [np.zeros_like(p) for p in self.params]
        self.v = [np.zeros_like(p) for p in self.params]
        self.t = 0

    def encode(self, X):
        H = np.tanh(X @ self.W1 + self.b1)
        return H @ self.W2 + self.b2

    def forward(self, X):
        H1 = np.tanh(X @ self.W1 + self.b1)
        Z = H1 @ self.W2 + self.b2
        H2 = np.tanh(Z @ self.W3 + self.b3)
        Y = H2 @ self.W4 + self.b4
        return H1, Z, H2, Y

    def step(self, X, lr=1e-3):
        n = len(X)
        H1, Z, H2, Y = self.forward(X)
        dY = 2.0 * (Y - X) / n
        gW4 = H2.T @ dY; gb4 = dY.sum(0)
        dH2 = (dY @ self.W4.T) * (1 - H2 ** 2)
        gW3 = Z.T @ dH2; gb3 = dH2.sum(0)
        dZ = dH2 @ self.W3.T
        gW2 = H1.T @ dZ; gb2 = dZ.sum(0)
        dH1 = (dZ @ self.W2.T) * (1 - H1 ** 2)
        gW1 = X.T @ dH1; gb1 = dH1.sum(0)
        grads = [gW1, gb1, gW2, gb2, gW3, gb3, gW4, gb4]
        self.t += 1
        b1, b2, eps = 0.9, 0.999, 1e-8
        for p, g, m, v in zip(self.params, grads, self.m, self.v):
            m[:] = b1 * m + (1 - b1) * g
            v[:] = b2 * v + (1 - b2) * g * g
            mh = m / (1 - b1 ** self.t)
            vh = v / (1 - b2 ** self.t)
            p -= lr * mh / (np.sqrt(vh) + eps)
        return float(((Y - X) ** 2).mean())

    def fit(self, X, epochs=200, batch=128, seed=0, lr=1e-3):
        r = np.random.default_rng(seed)
        loss = None
        for _ in range(epochs):
            idx = r.permutation(len(X))
            for k in range(0, len(X), batch):
                loss = self.step(X[idx[k:k + batch]], lr=lr)
        return loss


def run_aurora(seed: int, budget: int = BUDGET) -> Run:
    """MAP-Elites cuyo archivo NATIVO vive en el latente 2D del autoencoder.

    Fase 1: init aleatorio; entrenar AE con las curvas validas del init.
    Fase 2: MAP-Elites en el latente. Al 50% del presupuesto, re-entrenar con
    todo lo visto y RE-CODIFICAR el archivo (re-encoding de AURORA).
    Toda evaluacion se proyecta igualmente en el archivo de REFERENCIA, que es
    donde se comparan cobertura y QD-score entre configuraciones.
    """
    run = Run(seed, budget)
    arch = EliteArchive()
    sigs: list[np.ndarray] = []
    evs_kept: list[dict] = []
    ae = TinyAE(seed=seed)
    z_lo = z_hi = None

    def latent_cell(sig):
        z = ae.encode(sig[None, :])[0]
        i = int(np.clip((z[0] - z_lo[0]) / (z_hi[0] - z_lo[0] + 1e-12) * GRID,
                        0, GRID - 1))
        j = int(np.clip((z[1] - z_lo[1]) / (z_hi[1] - z_lo[1] + 1e-12) * GRID,
                        0, GRID - 1))
        return i, j

    # --- init ---
    while run.evals < min(INIT_EVALS, budget):
        topo, x = random_genome(run)
        ev = run.evaluate(topo, x)
        if ev and ev.get("valid"):
            sig = curve_signature(ev["P_eval"])
            sigs.append(sig)
            ev["_sig"] = sig
            evs_kept.append({"topology": ev["topology"], "x": ev["x"],
                             "rms_mm": ev["rms_mm"], "_sig": sig})
    X = np.asarray(sigs)
    ae.fit(X, epochs=150, seed=seed)
    Zl = ae.encode(X)
    z_lo = np.quantile(Zl, 0.01, axis=0)
    z_hi = np.quantile(Zl, 0.99, axis=0)
    for e in evs_kept:
        arch.offer(latent_cell(e["_sig"]), e)

    retrained = False
    while run.evals < budget:
        if not retrained and run.evals >= budget // 2:
            # re-entrenamiento AURORA con TODO lo almacenado + re-encoding
            X = np.asarray(sigs[-20000:])
            ae.fit(X, epochs=80, seed=seed + 1)
            Zl = ae.encode(X)
            z_lo = np.quantile(Zl, 0.01, axis=0)
            z_hi = np.quantile(Zl, 0.99, axis=0)
            old = list(arch.cells.values())
            arch = EliteArchive()
            for e in old:
                arch.offer(latent_cell(e["_sig"]), e)
            retrained = True
        for _ in range(BATCH):
            if run.evals >= budget:
                break
            if not arch.cells or run.rng.random() < IMMIGRANT_P:
                topo, x = random_genome(run)
            else:
                parent = arch.sample_parent(run.rng)
                topo, x = parent["topology"], mutate(run, topo := parent["topology"], parent["x"])
            ev = run.evaluate(topo, x)
            if ev and ev.get("valid"):
                sig = curve_signature(ev["P_eval"])
                sigs.append(sig)
                rec = {"topology": ev["topology"], "x": ev["x"],
                       "rms_mm": ev["rms_mm"], "_sig": sig}
                arch.offer(latent_cell(sig), rec)
    run.native_archive = arch
    run.native_space = "latent-AE"
    return run


# =====================================================================
# C — regla de reemplazo parametrizada
# =====================================================================
def make_replace_rule(theta: np.ndarray):
    """Regla logistica sobre [mejora_normalizada, vecinos_llenos, edad, misma
    familia]. theta[4] es el sesgo. La regla estandar equivale a aceptar solo
    con mejora > 0; esta puede aceptar retrocesos pequenos si, p.ej., la celda
    lleva mucho tiempo estancada. El mejor global se registra aparte, asi que
    esto no puede perder el mejor resultado."""
    th = np.asarray(theta, float)

    def rule(ev, inc, age, neighbors_filled):
        delta = (inc["rms_mm"] - ev["rms_mm"]) / max(inc["rms_mm"], 1e-9)
        feats = np.array([np.tanh(4.0 * delta),
                          neighbors_filled,
                          np.tanh(age / 40.0),
                          1.0 if ev["topology"] == inc["topology"] else 0.0])
        s = float(feats @ th[:4] + th[4])
        return 1.0 / (1.0 + np.exp(-s)) > 0.5
    return rule


def tune_replace_rule(n_theta=16, tune_seed=999, tune_budget=15000):
    """Busqueda aleatoria de theta en la semilla de ajuste. Se elige el mejor
    QD-score y se CONGELA antes de ver las semillas de evaluacion."""
    r = np.random.default_rng(tune_seed)
    results = []
    for k in range(n_theta):
        theta = np.concatenate([r.normal(0, 2.0, 4), r.normal(0, 1.0, 1)])
        run = run_map_elites(tune_seed, budget=tune_budget,
                             replace_rule=make_replace_rule(theta))
        m = run.metrics()
        results.append({"k": k, "theta": theta.tolist(),
                        "qd": m["qd_score"], "cov": m["coverage"],
                        "best": m["best_rms_mm"]})
        print(f"  theta[{k:02d}] qd={m['qd_score']:.2f} cov={m['coverage']:.3f}"
              f" best={m['best_rms_mm']:.3f}", flush=True)
    best = max(results, key=lambda d: d["qd"])
    return best, results


# =====================================================================
# E — pipeline existente instrumentado (sin modificar archivos)
# =====================================================================
def run_existing_pipeline(seed: int) -> dict:
    """Ejecuta fef.search_topology para las 8 familias con los valores por
    defecto del paper, contando CADA simulacion mediante un envoltorio del
    simulador en memoria (monkeypatch temporal; los archivos no se tocan).
    Toda evaluacion valida se proyecta al archivo de referencia."""
    run = Run(seed, budget=10 ** 9)          # presupuesto natural: sin tope
    orig_sim = fef.simulate_topology
    orig_eval = fef.evaluate_candidate
    sims = {"n": 0}

    def counting_sim(topo, x, theta):
        sims["n"] += 1
        return orig_sim(topo, x, theta)

    def logging_eval(topo, x, target, theta, allow_phase=True):
        ev = orig_eval(topo, x, target, theta, allow_phase=allow_phase)
        run.evals += 1
        if ev.get("valid"):
            d1, d2 = qd_core.descriptors(ev)
            run.ref.offer(ev, d1, d2)
            if run.best is None or ev["rms_mm"] < run.best["rms_mm"]:
                run.best = {k: ev[k] for k in
                            ("rms_mm", "max_mm", "topology", "accepted")}
                run.best["x"] = np.asarray(ev["x"], float).tolist()
            if (run.evals_to_pass is None and ev["rms_mm"] <= 1.0
                    and ev["max_mm"] <= 2.0):
                run.evals_to_pass = run.evals
        return ev

    fef.simulate_topology = counting_sim
    fef.evaluate_candidate = logging_eval
    try:
        for topo in FAMILIES:
            fef.search_topology(topo, run.target, run.theta)
    finally:
        fef.simulate_topology = orig_sim
        fef.evaluate_candidate = orig_eval

    m = run.metrics()
    m["sims"] = sims["n"]
    return m
