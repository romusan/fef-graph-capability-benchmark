"""Aplana las fuentes LaTeX para Editorial Manager.

EM no procesa subcarpetas: todos los ficheros deben quedar al mismo nivel. El
manuscrito usa `\\graphicspath{{figures/}}` e `\\input{tables/...}`, asi que hay
que copiar figuras y tablas a la raiz del paquete y reescribir las rutas.

No toca nada del proyecto: escribe una copia nueva en envio_RCIM/fuentes_em/.

Uso: py envio_RCIM/_aplana_para_em.py
"""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = Path(__file__).resolve().parent
PAPER = AQUI.parent / "paper"
DESTINO = AQUI / "fuentes_em"

TEX = PAPER / "FEF_Graph_Mechanism_Family_Q1.tex"


def main() -> None:
    if not TEX.exists():
        sys.exit(f"  *** no encuentro {TEX.name}")
    if DESTINO.exists():
        shutil.rmtree(DESTINO)
    DESTINO.mkdir(parents=True)

    s = TEX.read_text(encoding="utf-8")

    # Las tablas se incrustan: un \input menos es un fichero menos que perder.
    def mete_tabla(m):
        f = PAPER / "tables" / m.group(1)
        if not f.exists():
            sys.exit(f"  *** falta la tabla {m.group(1)}")
        return f.read_text(encoding="utf-8").rstrip()

    s, n_tab = re.subn(r"\\input\{tables/([^}]+)\}", mete_tabla, s)

    # Las figuras se copian a la raiz y se quita graphicspath.
    figs = set(re.findall(r"\\includegraphics\[[^\]]*\]\{([^}]+)\}", s))
    for f in sorted(figs):
        origen = PAPER / "figures" / f
        if not origen.exists():
            sys.exit(f"  *** falta la figura {f}")
        shutil.copy2(origen, DESTINO / f)
    s = re.sub(r"\\graphicspath\{\{figures/\}\}\s*\n", "", s)

    (DESTINO / TEX.name).write_text(s, encoding="utf-8")

    # Comprobacion: no puede quedar ninguna ruta con subcarpeta.
    restos = re.findall(r"\{(?:figures|tables)/[^}]*\}", s)
    if restos:
        sys.exit(f"  *** quedan rutas con subcarpeta: {restos[:3]}")

    print(f"  OK  {TEX.name}: {n_tab} tablas incrustadas, {len(figs)} figuras copiadas")
    print(f"  OK  sin rutas con subcarpeta")
    print(f"  -> {DESTINO}")
    print(f"     compilar ahi con pdflatex antes de subir, para confirmar que")
    print(f"     el aplanado no rompio nada")


if __name__ == "__main__":
    main()
