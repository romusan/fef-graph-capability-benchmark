"""Agrega P1 (8 familias x medium x 30 semillas) en la tabla del paper.

Por familia: pass-rate RMS<=1 con IC Wilson 95%, pass-rate estricto
(RMS<=1 y MAX<=2), mediana e IQR del mejor RMS por semilla, y tests exactos
de Fisher por pares frente a Watt I (con Holm dentro de esa familia de
comparaciones). Emite paper/tables/medium_all8.tex y un resumen por consola.

Uso:  py aggregate_p1.py     (cuando p1_log.txt diga "P1 COMPLETO")
"""
from __future__ import annotations

import glob
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BS = chr(92)
ROW = BS + BS

NICE = {"fourbar": "Four-bar", "slidercrank": "Slider-crank",
        "wattI": "Watt I", "wattII": "Watt II",
        "stephensonI": "Stephenson I", "stephensonII": "Stephenson II",
        "stephensonIII": "Stephenson III", "gearedFiveBar": "Geared five-bar"}
ORDEN = ["wattI", "stephensonIII", "stephensonI", "gearedFiveBar",
         "slidercrank", "fourbar", "wattII", "stephensonII"]


def wilson(k, n, z=1.959963984540054):
    if n == 0:
        return (0.0, 1.0)
    p = k / n
    d = 1 + z * z / n
    c = p + z * z / (2 * n)
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def holm(ps):
    orden = sorted(ps, key=ps.get)
    m, out, acc = len(orden), {}, 0.0
    for i, k in enumerate(orden):
        acc = max(acc, min(1.0, ps[k] * (m - i)))
        out[k] = acc
    return out


def main():
    filas = []
    for d in sorted(glob.glob(str(HERE / "p1_runs" / "*"))):
        arc = Path(d) / "candidate_archive.csv"
        if not arc.exists():
            continue
        seed = int(Path(d).name.split("seed")[1].split("_")[0])
        df = pd.read_csv(arc)
        for fam, g in df.groupby("topology"):
            # Misma regla que extract_family_rows() de run_paper_statistics.py,
            # el extractor que produjo las tablas por familia del paper.
            mejor = g.sort_values(["rms_mm", "max_mm"]).iloc[0]
            filas.append({"seed": seed, "family": fam,
                          "rms": float(mejor["rms_mm"]),
                          "mx": float(mejor["max_mm"])})
    df = pd.DataFrame(filas)
    n_seeds = df["seed"].nunique()
    print(f"semillas agregadas: {n_seeds}")

    res = []
    for fam in ORDEN:
        g = df[df.family == fam]
        n = len(g)
        k = int((g.rms <= 1.0).sum())
        ks = int(((g.rms <= 1.0) & (g.mx <= 2.0)).sum())
        lo, hi = wilson(k, n)
        res.append({"fam": fam, "n": n, "k": k, "ks": ks, "lo": lo, "hi": hi,
                    "med": g.rms.median(), "q1": g.rms.quantile(.25),
                    "q3": g.rms.quantile(.75)})

    # Fisher de cada familia frente a Watt I (pass RMS<=1), Holm
    kw = next(r for r in res if r["fam"] == "wattI")
    ps = {}
    for r in res:
        if r["fam"] == "wattI":
            continue
        tabla = [[kw["k"], kw["n"] - kw["k"]], [r["k"], r["n"] - r["k"]]]
        ps[r["fam"]] = fisher_exact(tabla)[1]
    adj = holm(ps)

    out = [r"\begin{table}[htbp]", r"\centering",
           r"\caption{Medium target, compact geometry, all eight families "
           r"under identical conditions and search configuration, "
           f"{n_seeds} seeds per family. Pass = RMS $" + BS + "leq$ 1 mm; "
           r"strict = RMS $\leq$ 1 mm and maximum error $\leq$ 2 mm. CI is "
           r"the Wilson 95\% interval for the pass rate. $p_{\mathrm{Holm}}$ "
           r"is the two-sided Fisher exact test of each family's pass count "
           r"against Watt I, Holm-corrected within this family of seven "
           r"comparisons.}",
           r"\label{tab:medium_all8}", r"\small",
           r"\begin{tabular}{lrrrrrl}", r"\hline",
           "Family & Pass & Strict & Pass rate CI & Median RMS & IQR & "
           "$p_{\\mathrm{Holm}}$ vs Watt I " + ROW, r"\hline"]
    for r in res:
        p = adj.get(r["fam"])
        ptxt = "---" if p is None else ("$<$0.001" if p < 1e-3 else f"{p:.3f}")
        out.append(
            f"{NICE[r['fam']]} & {r['k']}/{r['n']} & {r['ks']}/{r['n']} & "
            f"[{r['lo']:.2f}, {r['hi']:.2f}] & {r['med']:.3f} & "
            f"[{r['q1']:.3f}, {r['q3']:.3f}] & {ptxt} " + ROW)
        print(f"{NICE[r['fam']]:16s} pass {r['k']:2d}/{r['n']} "
              f"estricto {r['ks']:2d}/{r['n']} med {r['med']:7.3f} "
              f"CI[{r['lo']:.2f},{r['hi']:.2f}]"
              + ("" if p is None else f"  pHolm={p:.4f}"))
    out += [r"\hline", r"\end{tabular}", r"\end{table}"]
    dest = ROOT / "paper" / "tables" / "medium_all8.tex"
    dest.write_text("\n".join(out), encoding="utf-8")
    print(f"-> {dest}")


if __name__ == "__main__":
    main()
