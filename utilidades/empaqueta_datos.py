"""Reune las salidas por semilla de cada campana del manuscrito en datos/.

El apartado de disponibilidad del paper nombra siete campanas. Este script las
copia con nombres que dicen que experimento es cada una, en vez de la marca de
tiempo con que las escribio el motor, y genera un manifiesto con el recuento de
ficheros y una huella SHA-256 por campana.

No modifica nada del origen: solo lee.

Uso: py utilidades/empaqueta_datos.py <raiz-de-FEF-Graph>
"""
from __future__ import annotations

import hashlib
import json
import shutil
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")

# destino -> (origen relativo a la raiz, que tabla o seccion alimenta)
CAMPANAS = {
    "benchmark/easy_all_families": (
        "codigos/paper_stats_20260521_152557/runs",
        "Table 5, easy target, eight families, 30 seeds"),
    "benchmark/medium_wattI": (
        "codigos/paper_stats_20260521_121726/runs",
        "Tables 4 and 5, medium target under compact geometry, 30 seeds"),
    "benchmark/medium_hard_relaxed": (
        "codigos/paper_stats_20260521_131507/runs",
        "Table 7, medium-hard target under the relaxed 260 mm envelope"),
    "packaging/medium_hard_compact_wattI": (
        "codigos/c4_mh_compacto",
        "Table 7, same target under the compact 180 mm envelope, same seeds"),
    "packaging/medium_hard_geared_five_bar": (
        "codigos/c5_mh_familias",
        "Table 7, geared five-bar under both envelopes, 30 seeds per cell"),
    "policy/all_families_medium": (
        "qd_experiments/revision/p1_runs",
        "Table 9, benchmark policy applied to all eight families, 30 seeds"),
    "policy/qd_arms": (
        "qd_experiments/results",
        "Table 10, baseline, natural-budget and matched-budget arms, 15 seeds each"),
}


def huella(d: Path) -> tuple[str, int, int]:
    """SHA-256 sobre nombre y contenido de cada fichero, en orden estable."""
    h = hashlib.sha256()
    n = tam = 0
    for f in sorted(d.rglob("*")):
        if not f.is_file():
            continue
        h.update(str(f.relative_to(d)).replace("\\", "/").encode())
        b = f.read_bytes()
        h.update(b)
        n += 1
        tam += len(b)
    return h.hexdigest(), n, tam


def main() -> None:
    if len(sys.argv) < 2:
        sys.exit("  uso: py utilidades/empaqueta_datos.py <raiz-de-FEF-Graph>")
    raiz = Path(sys.argv[1]).resolve()
    destino = Path(__file__).resolve().parent.parent / "datos"
    if destino.exists():
        shutil.rmtree(destino)

    manifiesto = {"campanas": {}}
    total_f = total_b = 0
    for nombre, (origen, para_que) in CAMPANAS.items():
        src = raiz / origen
        if not src.is_dir():
            sys.exit(f"  *** no encuentro {src}")
        dst = destino / nombre
        shutil.copytree(src, dst)
        hh, n, tam = huella(dst)
        manifiesto["campanas"][nombre] = {
            "origen": origen, "alimenta": para_que,
            "ficheros": n, "bytes": tam, "sha256": hh}
        total_f += n
        total_b += tam
        print(f"  {nombre:40s} {n:5d} fich.  {tam/1024/1024:6.1f} MB  {hh[:12]}")

    manifiesto["total_ficheros"] = total_f
    manifiesto["total_bytes"] = total_b
    (destino / "MANIFIESTO.json").write_text(
        json.dumps(manifiesto, indent=2, ensure_ascii=False), encoding="utf-8")
    print(f"\n  total: {total_f} ficheros, {total_b/1024/1024:.1f} MB")
    print(f"  -> {destino}")


if __name__ == "__main__":
    main()
