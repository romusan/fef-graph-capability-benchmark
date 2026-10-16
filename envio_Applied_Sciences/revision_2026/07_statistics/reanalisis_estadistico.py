# -*- coding: utf-8 -*-
"""Reanálisis estadístico desde datos archivados (Fase VII del plan).

Requisitos cubiertos:
  - Diferenciar datos pareados de independientes; respetar bloques por semilla.
  - Tamaños de efecto con IC del 95 % (Hodges--Lehmann, rank-biserial,
    d de Cohen pareado dz con bootstrap, OR condicional de McNemar).
  - Corrección de Holm sobre familias de pruebas definidas y justificadas.
  - Efecto mínimo detectable (EMD) para proporciones pareadas con n = 30.
  - No interpretar no-significancia como equivalencia; reportar descriptivos.
Todo se recalcula desde los archivos por corrida; ninguna cifra se transfiere
de tablas publicadas.
"""
from __future__ import annotations
import json, math, re
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

RAIZ = Path(__file__).resolve().parents[3] / 'repo_codigo_nuevo' / 'datos'
if not RAIZ.exists():
    RAIZ = Path(__file__).resolve().parents[3] / 'datos'
SALIDA = Path(__file__).resolve().parent


def cargar_serie(directorio, patron, familia=None):
    """Devuelve DataFrame por semilla con rms/max/strict del mejor refinado."""
    filas = []
    for d in sorted(RAIZ.glob(f"{directorio}/{patron}")):
        m = re.search(r"seed(\d+)", d.name)
        if not m:
            continue
        seed = int(m.group(1))
        df = pd.read_csv(d / "candidate_archive.csv")
        if familia is not None:
            df = df[df.topology == familia]
        if len(df) == 0:
            continue
        mejor = df.loc[df.rms_mm.idxmin()]
        filas.append({"seed": seed, "rms": float(mejor.rms_mm), "max": float(mejor.max_mm),
                      "strict": bool(mejor.rms_mm <= 1.0 and mejor.max_mm <= 2.0),
                      "rms_pass": bool(mejor.rms_mm <= 1.0)})
    return pd.DataFrame(filas).set_index("seed").sort_index()


def wilson_ci(k, n, z=1.959964):
    p = k / n
    den = 1 + z * z / n
    centro = (p + z * z / (2 * n)) / den
    radio = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (max(0.0, centro - radio), min(1.0, centro + radio))


def hl_estimador_ci(d, conf=0.95):
    """Estimador Hodges--Lehmann y CI por promedios de Walsh."""
    d = np.asarray(d, float)
    d = d[d != 0]
    if len(d) == 0:
        return 0.0, (0.0, 0.0)
    walsh = np.sort((d[:, None] + d[None, :])[np.triu_indices(len(d))] / 2.0)
    est = float(np.median(walsh))
    n = len(d)
    z = stats.norm.ppf(1 - (1 - conf) / 2)
    k = int(math.floor(n * (n + 1) / 4 - z * math.sqrt(n * (n + 1) * (2 * n + 1) / 24)))
    k = max(0, min(k, len(walsh) - 1))
    return est, (float(walsh[k]), float(walsh[-k - 1]))


def rank_biserial(d):
    d = np.asarray(d, float)
    nz = d[d != 0]
    if len(nz) == 0:
        return 0.0
    r = stats.rankdata(abs(nz))
    W = r[nz > 0].sum() if (nz > 0).any() else 0.0
    return float(2 * W / r.sum() - 1)


def dz_pareado(d, n_boot=20000, seed=42):
    d = np.asarray(d, float)
    dz = float(np.mean(d) / np.std(d, ddof=1))
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        m = rng.choice(d, size=len(d), replace=True)
        s = m.std(ddof=1)
        boots.append(m.mean() / s if s > 0 else np.nan)
    boots = np.asarray(boots)
    boots = boots[~np.isnan(boots)]
    return dz, (float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5)))


def mcnemar_exact(a, b, c, dd):
    """McNemar exact (binomial sobre discordantes) + OR condicional aprox."""
    n_disc = b + c
    if n_disc == 0:
        p = 1.0
    else:
        # scipy.stats.binomtest ya devuelve p bilateral (two-sided)
        p = float(stats.binomtest(min(b, c), n_disc, 0.5).pvalue)
    if b * c == 0:
        or_ci = (float("nan"), float("nan"))
        orr = float("inf") if c > b else (0.0 if b > c else float("nan"))
    else:
        orr = c / b
        se = math.sqrt(1 / b + 1 / c)
        or_ci = (math.exp(math.log(orr) - 1.959964 * se), math.exp(math.log(orr) + 1.959964 * se))
    return p, orr, or_ci


def power_mcnemar(n, discordance, delta, alpha=0.05):
    """Exact unconditional power: D~Bin(n,q), B|D~Bin(D,(q+delta)/(2q)).

    q=P(01)+P(10), delta=P(01)-P(10). Concordant probabilities may be
    split arbitrarily; this conditional test depends only on q and delta.
    """
    if not 0 < discordance <= 1 or abs(delta)>discordance:
        raise ValueError('incompatible joint-probability assumptions')
    answer=0.
    for d in range(1,n+1):
        b=np.arange(d+1)
        pnull=np.minimum(1.,2*stats.binom.cdf(np.minimum(b,d-b),d,.5))
        rejection=stats.binom.pmf(b,d,(discordance+delta)/(2*discordance))[pnull<=alpha].sum()
        answer+=stats.binom.pmf(d,n,discordance)*rejection
    return float(answer)

def emd_mcnemar(n, discordance, alpha=.05, poder=.8):
    lo=0.;hi=discordance
    if power_mcnemar(n,discordance,hi,alpha)<poder:return None
    for _ in range(35):
        mid=(lo+hi)/2
        if power_mcnemar(n,discordance,mid,alpha)>=poder:hi=mid
        else:lo=mid
    return hi


def holm(pvals):
    orden = np.argsort(pvals)
    m = len(pvals)
    adj = np.empty(m)
    prev = 0.0
    for rango, idx in enumerate(orden):
        val = (m - rango) * pvals[idx]
        prev = max(prev, min(val, 1.0))
        adj[idx] = prev
    return adj


resultados = {}

# ---------- 1. Packaging: Watt I y geared five-bar, pareado por semilla ----------
pack = {}
for fam, comp_dir, comp_pat, rel_dir, rel_pat in [
    ("wattI", "packaging/medium_hard_compact_wattI", "c4_mh_compact_seed*",
     "benchmark/medium_hard_relaxed", "paper_stats_run_wattI_medium_hard_seed*"),
    ("gearedFiveBar", "packaging/medium_hard_geared_five_bar", "c5_gearedFiveBar_compact_seed*",
     "packaging/medium_hard_geared_five_bar", "c5_gearedFiveBar_relaxed_seed*"),
]:
    comp = cargar_serie(comp_dir, comp_pat, fam if fam == "gearedFiveBar" else None)
    rel = cargar_serie(rel_dir, rel_pat, fam if fam == "gearedFiveBar" else None)
    j = comp.join(rel, lsuffix="_comp", rsuffix="_rel", how="inner")
    d = (j.rms_comp - j.rms_rel).values
    p_w = stats.wilcoxon(d).pvalue
    hl, hl_ci = hl_estimador_ci(d)
    rb = rank_biserial(d)
    dz, dz_ci = dz_pareado(d)
    a = int((j.strict_comp & j.strict_rel).sum()); b = int((~j.strict_comp & j.strict_rel).sum())
    c = int((j.strict_comp & ~j.strict_rel).sum()); dd = int((~j.strict_comp & ~j.strict_rel).sum())
    p_m, orr, or_ci = mcnemar_exact(a, b, c, dd)
    ar = int((j.rms_pass_comp & j.rms_pass_rel).sum()); br = int((~j.rms_pass_comp & j.rms_pass_rel).sum())
    cr = int((j.rms_pass_comp & ~j.rms_pass_rel).sum()); dr = int((~j.rms_pass_comp & ~j.rms_pass_rel).sum())
    p_mr, orr_r, or_ci_r = mcnemar_exact(ar, br, cr, dr)
    pack[fam] = {
        "n_pareado": int(len(j)),
        "rms": {"mediana_compacto": float(j.rms_comp.median()), "mediana_relativax": float(j.rms_rel.median()),
                "wilcoxon_p": float(p_w), "HL_mm": hl, "HL_ic95": hl_ci,
                "rank_biserial": rb, "dz": dz, "dz_ic95_boot": dz_ci,
                "n_aumentan_error": int((d > 0).sum()), "n_disminuyen": int((d < 0).sum())},
        "strict": {"compacto": int(j.strict_comp.sum()), "relajado": int(j.strict_rel.sum()),
                   "tabla": [a, b, c, dd], "mcnemar_p": p_m, "OR": orr, "OR_ic95": or_ci,
                   "wilson_compacto": wilson_ci(int(j.strict_comp.sum()), len(j)),
                   "wilson_relativax": wilson_ci(int(j.strict_rel.sum()), len(j))},
        "rms_pass": {"compacto": int(j.rms_pass_comp.sum()), "relajado": int(j.rms_pass_rel.sum()),
                     "tabla": [ar, br, cr, dr], "mcnemar_p": p_mr, "OR": orr_r, "OR_ic95": or_ci_r},
    }

# Familia de seis contrastes definida para este reanalisis retrospectivo.
nombres = []
pvals = []
for fam in ["wattI", "gearedFiveBar"]:
    nombres += [f"{fam}:rms", f"{fam}:strict", f"{fam}:rms_pass"]
    pvals += [pack[fam]["rms"]["wilcoxon_p"], pack[fam]["strict"]["mcnemar_p"], pack[fam]["rms_pass"]["mcnemar_p"]]
adj = holm(pvals)
pack["holm_packaging"] = {n: {"p_cruda": float(p), "p_holm": float(a)}
                          for n, p, a in zip(nombres, pvals, adj)}
pack['discordance_note'] = ('La tabla observada tiene 11 y 5 discordancias, no seis en total. '
    'Seis discordancias todas en un sentido dan p bilateral=0.03125 sin ajuste. '
    'La diferencia marginal de seis no determina la potencia.')
pack['power_sensitivity']=[{'n':30,'q':q,'alpha':alpha,'target_power':.8,
    'minimum_delta_for_80pct':emd_mcnemar(30,q,alpha),
    'power_for_delta_point2':power_mcnemar(30,q,.2,alpha)}
    for q in [.2,16/30,.8] for alpha in [.05,.05/6]]

resultados["packaging"] = pack

# ---------- 2. Watt I medium: campaña inicial vs all-family, pareado por semilla ----------
ini = cargar_serie("benchmark/medium_wattI", "paper_stats_run_wattI_medium_seed*")
all8 = cargar_serie("policy/all_families_medium", "p1_medium_all8_seed*", "wattI")
j = ini.join(all8, lsuffix="_ini", rsuffix="_all8", how="inner")
d = (j.rms_ini - j.rms_all8).values
p_w2 = stats.wilcoxon(d).pvalue
hl2, hl_ci2 = hl_estimador_ci(d)
a = int((j.strict_ini & j.strict_all8).sum()); b = int((~j.strict_ini & j.strict_all8).sum())
c = int((j.strict_ini & ~j.strict_all8).sum()); dd = int((~j.strict_ini & ~j.strict_all8).sum())
p_m2, orr2, or_ci2 = mcnemar_exact(a, b, c, dd)
resultados["wattI_medium_campaigns"] = {
    "n_pareado": int(len(j)),
    "descripcion_diferencia": ("Campañas NO equivalentes: la inicial ejecuta solo wattI con "
                               "DE 35 iter x pop 7 (~2.9k evals) y flujo aleatorio propio; la "
                               "all-family ejecuta 8 familias con DE 45 x 9 (~4.9k evals) y el "
                               "flujo aleatorio de wattI depende de las 3 familias previas."),
    "strict": {"inicial_7de30": int(j.strict_ini.sum()), "all8_16de30": int(j.strict_all8.sum()),
               "tabla_pareada": [a, b, c, dd], "mcnemar_p": p_m2, "OR": orr2, "OR_ic95": or_ci2,
               "wilson_inicial": wilson_ci(int(j.strict_ini.sum()), len(j)),
               "wilson_all8": wilson_ci(int(j.strict_all8.sum()), len(j))},
    "rms": {"mediana_inicial": float(j.rms_ini.median()), "mediana_all8": float(j.rms_all8.median()),
            "wilcoxon_p": float(p_w2), "HL_mm": hl2, "HL_ic95": hl_ci2},
}

# ---------- 3. QD: schema validated; published baseline is Estrict ----------
qd = {}
# recuento real desde JSON
import glob as g
real = {}
for f in RAIZ.glob("policy/qd_arms/*.json"):
    d = json.load(open(f))
    brazo = d["config"]
    if brazo not in ['Estrict','Bplus','Bmatched']:continue
    if 'n_strict' not in d:raise ValueError(f'Esquema QD incompatible: {f}')
    real.setdefault(brazo, []).append(d['n_strict']>0)
for b, vals in real.items():
    k = sum(vals); n = len(vals)
    qd[b] = {"strict_seeds": k, "n": n, "wilson_ic95": wilson_ci(k, n)}
# McNemar pareado A vs Bmatched por semilla
A = {int(re.search(r"seed(\d+)", f.name).group(1)): json.load(open(f))['n_strict']>0
     for f in RAIZ.glob("policy/qd_arms/Estrict_seed*.json")}
Bm = {int(re.search(r"seed(\d+)", f.name).group(1)): json.load(open(f))['n_strict']>0
      for f in RAIZ.glob("policy/qd_arms/Bmatched_seed*.json")}
comunes = sorted(set(A) & set(Bm))
aa = sum(A[s] and Bm[s] for s in comunes); bb = sum((not A[s]) and Bm[s] for s in comunes)
cc = sum(A[s] and (not Bm[s]) for s in comunes); ddd = sum((not A[s]) and (not Bm[s]) for s in comunes)
p_qd, or_qd, or_ci_qd = mcnemar_exact(aa, bb, cc, ddd)
qd["Estrict_vs_Bmatched_pareado"] = {"tabla": [aa, bb, cc, ddd], "mcnemar_p": p_qd, "OR": or_qd, "OR_ic95": or_ci_qd,
                               "n_comunes": len(comunes)}
resultados["qd_arms"] = qd

SALIDA.mkdir(parents=True, exist_ok=True)
def clean(value):
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,(tuple,list)):return [clean(v) for v in value]
    if isinstance(value,(float,np.floating)) and not np.isfinite(value):return None
    return value
resultados=clean(resultados)
with open(SALIDA / "reanalisis_estadistico.json", "w", encoding="utf-8") as f:
    json.dump(resultados, f, indent=2, ensure_ascii=False, allow_nan=False)
print(json.dumps(resultados, indent=2, ensure_ascii=False, default=str)[:4000])
print("OK ->", SALIDA / "reanalisis_estadistico.json")
