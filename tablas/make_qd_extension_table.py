"""Genera paper/tables/qd_extension.tex desde los JSON de las campanas
Estrict (linea base) y Bplus (extension QD), 15 semillas por brazo.

Con 15 semillas la tabla por semilla tendria 30 filas y ademas enterraria el
resultado: las tasas de exito son indistinguibles y lo que difiere es que
familia encuentra cada politica. La tabla resume por desenlace y deja el
detalle por semilla a los datos publicados.
"""
import json
from pathlib import Path
from math import sqrt
from scipy.stats import fisher_exact

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BS = chr(92)
FIN = BS + BS
SEMILLAS = [101, 202, 303, 404, 505, 606, 707, 808, 909,
            1010, 1111, 1212, 1313, 1414, 1515]
N = len(SEMILLAS)


def carga(cfg):
    return {s: json.loads((HERE / "results" / f"{cfg}_seed{s}.json")
                          .read_text(encoding="utf-8")) for s in SEMILLAS}


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


E, B, M = carga("Estrict"), carga("Bplus"), carga("Bmatched")

CASOS = [("Seeds with at least one strict design", None),
         ("~~~from geared five-bar", "gearedFiveBar"),
         ("~~~from Watt I", "wattI"),
         ("~~~from Stephenson III", "stephensonIII")]

filas, crudos = [], []
for nombre, fam in CASOS:
    def cuenta(d, fam=fam):
        if fam is None:
            return sum(1 for s in SEMILLAS if d[s].get("n_strict", 0) > 0)
        return sum(1 for s in SEMILLAS
                   if fam in (d[s].get("strict_families") or []))
    a, b, m = cuenta(E), cuenta(B), cuenta(M)
    # El contraste que importa es el de presupuesto igualado contra la
    # linea base: es el unico controlado por coste. Son cuatro hipotesis
    # sobre los mismos datos, asi que se reporta tambien el p ajustado por
    # Holm; sin el, un 0.035 se leeria como significativo cuando no lo es.
    p = fisher_exact([[m, N - m], [a, N - a]])[1]
    lm, hm = wilson(m, N)
    crudos.append(p)
    filas.append([nombre, a, b, m, lm, hm, p])

# Holm sobre las cuatro comparaciones de esta tabla.
orden = sorted(range(len(crudos)), key=lambda i: crudos[i])
prev, ajust = 0.0, [0.0] * len(crudos)
for k, i in enumerate(orden):
    prev = max(prev, min(1.0, (len(crudos) - k) * crudos[i]))
    ajust[i] = prev
fmt = lambda x: "$<$0.001" if x < 0.001 else f"{x:.3f}"
filas = [f"{r[0]} & {r[1]}/{N} & {r[2]}/{N} & {r[3]}/{N} & "
         f"[{r[4]:.2f}, {r[5]:.2f}] & {fmt(r[6])} & {fmt(ajust[i])} " + FIN
         for i, r in enumerate(filas)]

tot = lambda d, k: sum(d[s][k] for s in SEMILLAS)
mejor = lambda d: min(d[s]["best_rms_mm"] for s in SEMILLAS)
se, sb, sm = (tot(d, "sims_total") for d in (E, B, M))
resumen = [
    f"Strict designs, total & {tot(E,'n_strict')} & {tot(B,'n_strict')} & "
    f"{tot(M,'n_strict')} & & --- & --- " + FIN,
    f"Simulations, total & {se:,} & {sb:,} & {sm:,} & & --- & --- " + FIN,
    f"Best RMS of any seed (mm) & {mejor(E):.3f} & {mejor(B):.3f} & "
    f"{mejor(M):.3f} & & --- & --- " + FIN,
]

out = [
    r"\begin{table}[htbp]", r"\centering",
    r"\caption{Search policies on the medium target under compact "
    r"geometry, " + str(N) + r" seeds each, same simulator, same acceptance "
    r"criteria and same least-squares refinement stage. The baseline draws "
    r"120 random samples and 16 shape-seeded templates per family with "
    r"differential evolution disabled; the extension couples a MAP-Elites "
    r"archive with surrogate-guided variation. Strict criterion: RMS $\leq$ "
    r"1 mm and maximum error $\leq$ 2 mm. Intervals are Wilson 95\% intervals "
    r"for the seed-wise rate of the approximately budget-matched extension. The raw $p$ is "
    r"the two-sided Fisher exact test of that arm against the baseline; the "
    r"Holm column adjusts the four comparisons of this table for multiplicity, "
    r"and none of them remains below 0.05 after adjustment. The "
    r"extension is reported twice: at its natural budget and at a search "
    r"budget reduced from 60{,}000 to 34{,}000 evaluations so that total "
    r"simulations come within 7\% of the baseline. All simulation counts are "
    r"measured.}",
    r"\label{tab:qd_extension}", r"\scriptsize",
    r"\resizebox{0.95\textwidth}{!}{%",
    r"\begin{tabular}{lrrrlrr}", r"\hline",
    r"& & \multicolumn{2}{c}{QD + surrogate} & & \multicolumn{2}{c}{$p$} " + FIN,
    r"\cline{3-4} \cline{6-7}",
    r"Outcome & Baseline & natural & matched & 95\% CI (matched) & raw & Holm " + FIN,
    r"\hline",
] + filas + [r"\hline"] + resumen + [
    r"\hline", r"\end{tabular}}", r"\end{table}", ""]

(ROOT / "paper" / "tables" / "qd_extension.tex").write_text(
    "\n".join(out), encoding="utf-8")
print(f"tabla generada con {N} semillas por brazo")
for f in filas:
    print("  " + f.replace(FIN, "").replace("~~~", "  "))
