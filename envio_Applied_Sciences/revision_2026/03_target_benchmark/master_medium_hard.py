# -*- coding: utf-8 -*-
"""Maestro de la campana medium-hard x 8 familias (etapa 2 del PLAN_EXPERIMENTAL).

Cola: 8 familias x 30 semillas x perfiles (compact + relaxed) = 480 corridas.
Orden: todas las compact primero, luego las relaxed (el perfil compact, que
alimenta la respuesta a los revisores, termina antes).

Paralelismo: 3 procesos worker simultaneos (maquina de 4 nucleos logicos /
16 GB RAM). Cada worker escribe su propia carpeta; el maestro es reanudable:
al arrancar filtra las corridas cuyo run_summary.json ya existe.

Un worker que falla (exit code != 0 o sin run_summary.json) se reintenta
hasta 2 veces; si sigue fallando se registra en failures.json y la campana
continua con las demas.

Progreso en progreso.json: hechas/pendientes/fallidas + tiempos medios y ETA.
El .bat corre_campana.bat reinicia este maestro automaticamente si se cae;
para detenerlo de forma limpia crear el archivo STOP_CAMPANA.txt en esta
carpeta.

Uso:
  python master_medium_hard.py                # compact + relaxed (480)
  python master_medium_hard.py --profiles compact   # solo compact (240)
  python master_medium_hard.py --processes 2  # cambiar paralelismo
"""
from __future__ import annotations
import argparse, json, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

AQUI = Path(__file__).resolve().parent
WORKER = AQUI / "worker_medium_hard.py"
FAMILIAS = ["fourbar", "slidercrank", "wattI", "wattII",
            "stephensonI", "stephensonII", "stephensonIII", "gearedFiveBar"]
SEEDS = range(1, 31)
MAX_REINTENTOS = 2


def carpeta(perfil: str, familia: str, semilla: int) -> Path:
    return AQUI / f"medium_hard_{perfil}" / f"{familia}_seed{semilla:03d}"


def construir_cola(perfiles: list[str]) -> list[tuple[str, str, int]]:
    cola = []
    for perfil in perfiles:
        for familia in FAMILIAS:
            for semilla in SEEDS:
                cola.append((perfil, familia, semilla))
    return cola


def correr(tarea: tuple[str, str, int]) -> tuple[tuple[str, str, int], bool, str]:
    perfil, familia, semilla = tarea
    for intento in range(1, MAX_REINTENTOS + 2):  # 1 intento + 2 reintentos
        proc = subprocess.run(
            [sys.executable, str(WORKER), familia, str(semilla), perfil],
            capture_output=True, text=True, timeout=7200, cwd=str(AQUI))
        salida_ok = proc.returncode == 0 and (carpeta(perfil, familia, semilla) / "run_summary.json").exists()
        if salida_ok:
            linea = (proc.stdout or "").strip().splitlines()
            return tarea, True, linea[-1] if linea else "ok"
        if intento <= MAX_REINTENTOS:
            time.sleep(5)
    return tarea, False, (proc.stderr or proc.stdout or "sin salida")[-400:]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--profiles", nargs="+", default=["compact", "relaxed"],
                    choices=["compact", "relaxed"])
    ap.add_argument("--processes", type=int, default=3)
    args = ap.parse_args()

    cola_total = construir_cola(args.profiles)
    pendientes = [t for t in cola_total
                  if not (carpeta(*t) / "run_summary.json").exists()]
    hechas = len(cola_total) - len(pendientes)

    print(f"[maestro] perfiles={args.profiles} total={len(cola_total)} "
          f"hechas={hechas} pendientes={len(pendientes)} procesos={args.processes}")
    if not pendientes:
        print("[maestro] nada pendiente, campana completa.")
        return 0

    t_campana = time.time()
    inicio_sesion = time.strftime("%Y-%m-%d %H:%M:%S")
    tiempos: list[float] = []
    fallidas: list[dict] = []
    completadas = 0

    with ThreadPoolExecutor(max_workers=args.processes) as ex:
        futuros = {ex.submit(correr, t): t for t in pendientes}
        for fut in as_completed(futuros):
            tarea, ok, detalle = fut.result()
            perfil, familia, semilla = tarea
            completadas += 1
            if ok:
                # tiempo de la corrida desde el resumen del worker
                try:
                    r = json.loads((carpeta(*tarea) / "run_summary.json").read_text(encoding="utf-8"))
                    tiempos.append(float(r.get("elapsed_seconds", 0.0)))
                except Exception:
                    pass
                print(f"[{completadas}/{len(pendientes)}] OK {perfil}/{familia} s{semilla}: {detalle}")
            else:
                fallidas.append({"perfil": perfil, "familia": familia, "semilla": semilla,
                                 "detalle": detalle})
                print(f"[{completadas}/{len(pendientes)}] FALLO {perfil}/{familia} s{semilla} "
                      f"tras {MAX_REINTENTOS} reintentos", file=sys.stderr)
            # progreso persistente
            t_medio = sum(tiempos) / len(tiempos) if tiempos else 580.0
            restantes = len(pendientes) - completadas
            eta_s = restantes * t_medio / args.processes
            progreso = {
                "perfiles": args.profiles, "total_cola": len(cola_total),
                "hechas_previas": hechas, "completadas_esta_sesion": completadas,
                "pendientes": restantes, "fallidas": len(fallidas),
                "tiempo_medio_s": round(t_medio, 1),
                "eta_minutos": round(eta_s / 60, 1),
                "sesion_inicio": inicio_sesion,
                "sesion_tiempo_transcurrido_min": round((time.time() - t_campana) / 60, 1),
            }
            (AQUI / "progreso.json").write_text(json.dumps(progreso, indent=2), encoding="utf-8")

    if fallidas:
        (AQUI / "failures.json").write_text(json.dumps(fallidas, indent=2), encoding="utf-8")
        print(f"[maestro] terminado con {len(fallidas)} fallidas -> failures.json")
        return 1
    print("[maestro] campana completa, sin fallidas.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
