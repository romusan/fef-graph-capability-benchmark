"""Re-simula en proceso independiente el Watt II estricto de la Tabla 9.

La Tabla 9 reporta Watt II con 1/30 estricto en el objetivo medio bajo
geometria compacta. Si ese diseno reproduce sus errores, el numero de
familias con al menos un diseno estricto pasa de tres a cuatro y la
conclusion de 5.6 debe corregirse.

El vector sale de candidate_archive.csv de la semilla 026. Orden de variables
segun la Tabla 1 del manuscrito para Watt II: L0, L2, L3, L4, r_c, phi_c,
gx, gy, L5, L6, pu, pv (12 variables). Se evalua con el mismo simulador y el
mismo objetivo medium, sin refinar.

Uso: py qd_experiments/revision/verifica_wattII.py
"""
from __future__ import annotations

import csv
import glob
import json
import os
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
AQUI = Path(__file__).resolve().parent
sys.path.insert(0, str(AQUI.parent))

os.environ["FEF_GEOMETRY_PROFILE"] = "compact"
os.environ["FEF_MAX_PRINTABLE_LINK_MM"] = "180.0"
os.environ["FEF_PHASE_STRIDE"] = "1"

import numpy as np
import qd_core
import fef_mechanism_family_synthesis as fef

ORDEN_WATTII = ["L0", "L2", "L3", "L4", "r_c", "phi_c",
                "gx", "gy", "L5", "L6", "pu", "pv"]

d = glob.glob(str(AQUI / "p1_runs" / "p1_medium_all8_seed026_*"))[0]
with open(os.path.join(d, "candidate_archive.csv"), encoding="utf-8") as f:
    filas = [r for r in csv.DictReader(f) if r["topology"] == "wattII"]
mejor = min(filas, key=lambda r: (float(r["rms_mm"]), float(r["max_mm"])))
arch_rms, arch_max = float(mejor["rms_mm"]), float(mejor["max_mm"])

faltan = [k for k in ORDEN_WATTII if not mejor.get(k, "").strip()]
if faltan:
    sys.exit(f"  *** faltan variables en el archivo: {faltan}")
x = np.array([float(mejor[k]) for k in ORDEN_WATTII], dtype=float)

print(f"  semilla 026, Watt II, {len(x)} variables")
print(f"  archivado  : RMS {arch_rms:.4f} mm   MAX {arch_max:.4f} mm")

target, theta, meta = qd_core.target_and_theta()
ev = fef.evaluate_candidate("wattII", x, target, theta)
print(f"  valido     : {ev.get('valid')}")
if ev.get("valid"):
    r, m = float(ev["rms_mm"]), float(ev["max_mm"])
    print(f"  re-simulado: RMS {r:.4f} mm   MAX {m:.4f} mm")
    print(f"  coincide   : RMS {abs(r-arch_rms) < 5e-3}   MAX {abs(m-arch_max) < 5e-3}")
    print(f"  estricto   : {r <= 1.0 and m <= 2.0}")
    salida = {"topology": "wattII", "seed": 26, "variables": ORDEN_WATTII,
              "x": x.tolist(), "archivado": {"rms_mm": arch_rms, "max_mm": arch_max},
              "resimulado": {"rms_mm": r, "max_mm": m},
              "estricto": bool(r <= 1.0 and m <= 2.0)}
else:
    salida = {"topology": "wattII", "seed": 26, "valido": False,
              "motivo": {k: ev.get(k) for k in ("reason", "why", "status") if k in ev}}
    print(f"  *** el candidato no es valido al re-simular: {salida['motivo']}")

(AQUI / "wattII_verificacion.json").write_text(
    json.dumps(salida, indent=2, default=str), encoding="utf-8")
print(f"  -> wattII_verificacion.json")
