"""Genera paper/tables/packaging_effect.tex: diseno 2x2 del medium-hard.

Dos familias (Watt I, geared five-bar) por dos envolventes (relajado 260 mm,
compacto 180 mm), 30 semillas por celda. Las cuatro campanas comparten
objetivo, presupuesto de muestreo, evolucion diferencial y refinado; sus
run_summary.json difieren solo en topologia y en los dos campos del
envolvente.

Uso: py codigos/make_packaging_table.py
"""
from __future__ import annotations

import glob
import json
import os
import statistics as st
from pathlib import Path

from scipy.stats import binomtest, fisher_exact, wilcoxon

AQUI = Path(__file__).resolve().parent
ROOT = AQUI.parent
BS = chr(92)
FIN = BS + BS

FUENTES = {
    ("wattI", "relaxed"): ("paper_stats_20260521_131507/runs/"
                           "paper_stats_run_wattI_medium_hard_seed*/best_accuracy.json"),
    ("wattI", "compact"): "c4_mh_compacto/*/best_accuracy.json",
    ("g5b", "relaxed"): "c5_mh_familias/*_relaxed_seed*/best_accuracy.json",
    ("g5b", "compact"): "c5_mh_familias/*_compact_seed*/best_accuracy.json",
}
NOMBRE = {"wattI": "Watt I", "g5b": "Geared five-bar"}


def carga(patron):
    """Devuelve {semilla: (rms, max)}. Las cuatro campanas usan las
    semillas 1--30, asi que las celdas estan PAREADAS y los contrastes entre
    envolventes deben ser pareados: McNemar exacto para las tasas y Wilcoxon
    de rangos con signo para el error."""
    out = {}
    for p in sorted(glob.glob(str(AQUI / patron))):
        s = int(os.path.basename(os.path.dirname(p)).split("seed")[1][:3])
        m = json.load(open(p, encoding="utf-8")).get("metrics", {})
        if "rms_mm" in m:
            out[s] = (m["rms_mm"], m["max_mm"])
    return out


D = {k: carga(v) for k, v in FUENTES.items()}
n = len(D[("wattI", "relaxed")])

def mcnemar(r, c, ss, crit):
    """McNemar exacto sobre los pares discordantes."""
    b = sum(1 for s in ss if crit(r[s]) and not crit(c[s]))
    d = sum(1 for s in ss if crit(c[s]) and not crit(r[s]))
    return (binomtest(b, b + d, 0.5).pvalue if b + d else 1.0), b, d


filas = []
for fam in ("wattI", "g5b"):
    r, c = D[(fam, "relaxed")], D[(fam, "compact")]
    ss = sorted(set(r) & set(c))
    filas.append(BS + "textit{" + NOMBRE[fam] + "} & & & " + FIN)
    for etq, crit in [("~~~Passes RMS $" + BS + "leq$ 1 mm",
                       lambda t: t[0] <= 1.0),
                      ("~~~Meets strict criterion",
                       lambda t: t[0] <= 1.0 and t[1] <= 2.0)]:
        a = sum(1 for s in ss if crit(r[s]))
        b = sum(1 for s in ss if crit(c[s]))
        p, disc_b, disc_d = mcnemar(r, c, ss, crit)
        filas.append(f"{etq} & {a}/{len(ss)} & {b}/{len(ss)} & {p:.2f} " + FIN)
    pw = wilcoxon([r[s][0] for s in ss], [c[s][0] for s in ss]).pvalue
    ptxt = "$<$0.001" if pw < 0.001 else f"{pw:.3f}"
    filas.append(f"~~~Median RMS error (mm) & {st.median(r[s][0] for s in ss):.3f} & "
                 f"{st.median(c[s][0] for s in ss):.3f} & {ptxt} " + FIN)

# Contraste ENTRE familias dentro de cada envolvente: aqui las muestras son
# independientes (familias distintas), asi que Fisher es la prueba correcta.
entre = []
for perfil in ("relaxed", "compact"):
    a = sum(1 for x, y in D[("wattI", perfil)].values() if x <= 1.0 and y <= 2.0)
    b = sum(1 for x, y in D[("g5b", perfil)].values() if x <= 1.0 and y <= 2.0)
    p = fisher_exact([[a, n - a], [b, n - b]])[1]
    entre.append(f"{perfil}: {p:.1e}")

out = [
    r"\begin{table}[htbp]", r"\centering",
    r"\caption{Medium-hard target under two packaging envelopes and two "
    r"families, " + str(n) + r" seeds per cell. The four campaigns share the "
    r"target curve, the sampling budget, the differential-evolution stage and "
    r"the refinement depth; their recorded configurations differ only in the "
    r"family and in the two fields that define the envelope. Strict criterion: "
    r"RMS $\leq$ 1 mm and maximum error $\leq$ 2 mm. The two envelopes were run "
    r"on the same seeds 1--30, so the comparison within a family is paired: $p$ "
    r"is the exact McNemar test for the rates and the Wilcoxon signed-rank test "
    r"for the error. Between families within an envelope the samples are "
    r"independent and "
    r"the strict rates differ, by Fisher exact test, at $p=" + entre[0].split(": ")[1].replace("e-", r"\times 10^{-")
    + r"}$ (relaxed) and $p=" + entre[1].split(": ")[1].replace("e-", r"\times 10^{-")
    + r"}$ (compact).}",
    r"\label{tab:packaging_effect}", r"\small",
    r"\begin{tabular}{lrrr}", r"\hline",
    r"Outcome & Relaxed, 260 mm & Compact, 180 mm & $p$ " + FIN,
    r"\hline",
] + filas + [r"\hline", r"\end{tabular}", r"\end{table}", ""]

txt = "\n".join(out).replace(r"\times 10^{-0", r"\times 10^{-")
(ROOT / "paper" / "tables" / "packaging_effect.tex").write_text(txt, encoding="utf-8")
print(f"tabla 2x2 generada, {n} semillas por celda")
for f in filas:
    print("  " + f.replace(FIN, "").replace("~~~", "   "))
print("  entre familias: " + "  ".join(entre))
