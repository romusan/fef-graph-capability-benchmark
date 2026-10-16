"""
FEF multi-topology inverse synthesis + additive-manufacturing prototype exporter
================================================================================

Goal
----
Search three planar mechanism families for a moderate-accuracy desktop path-generation
application, then export editable CAD files (STEP) and STL files for the most printable/best candidate.

Families:
  1) fourbar
  2) slider-crank
  3) Stephenson III, using the exact two-loop geometric convention:
       Loop 1: O2-A-B-O4
       Loop 2: C-D-O4
     with C fixed in the local AB frame and output point P in the local CD frame.

Application framing
-------------------
This script intentionally avoids a high-precision welding claim. The default application is a
low-cost additive-manufactured desktop curve-tracing / inspection / dispensing demonstrator.
The default tolerance is therefore moderate:
    RMS <= 1.5 mm and MAX <= 3.0 mm

Outputs
-------
A timestamped folder containing:
  - target_path.csv
  - candidate_archive.csv
  - best_accuracy.json
  - best_printable.json
  - figures/*.png
  - stl_best_printable/*.stl
  - step_best_printable/*.step
  - BOM_and_assembly_notes.txt
  - results_package.zip

Dependencies
------------
Required: numpy, scipy, pandas, matplotlib
Optional but recommended: numba, shapely, trimesh, cadquery

Install:
    python -m pip install numpy scipy pandas matplotlib numba shapely trimesh

For editable STEP export in SolidWorks, install CadQuery:
    python -m pip install cadquery

Notes on STL/STEP generation
----------------------------
The STL exporter uses shapely for 2D geometry and a custom polygon extrusion routine,
so it does NOT require boolean operations or external CAD kernels.
The STEP exporter uses CadQuery/OpenCascade when available. STEP files are preferred
for SolidWorks editing because they preserve CAD B-rep solids better than STL meshes.
"""
from __future__ import annotations

import json
import math
import os
import platform
import time
import zipfile
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import differential_evolution, minimize, least_squares

try:
    import fef_cpp_core as cpp_core
    CPP_CORE_AVAILABLE = bool(cpp_core.CPP_CORE_AVAILABLE)
    CPP_CORE_IMPORT_ERROR = cpp_core.IMPORT_ERROR
except Exception as exc:
    cpp_core = None
    CPP_CORE_AVAILABLE = False
    CPP_CORE_IMPORT_ERROR = f"{type(exc).__name__}: {exc}"

try:
    from scipy.stats import qmc
    QMC_AVAILABLE = True
except Exception:
    qmc = None
    QMC_AVAILABLE = False

try:
    from numba import njit
    NUMBA_AVAILABLE = True
except Exception:
    NUMBA_AVAILABLE = False
    def njit(*args, **kwargs):
        def wrapper(f): return f
        return wrapper

try:
    from shapely.geometry import Point, LineString, Polygon, MultiPoint, box
    from shapely.ops import unary_union, triangulate
    SHAPELY_AVAILABLE = True
except Exception:
    SHAPELY_AVAILABLE = False

try:
    import trimesh
    TRIMESH_AVAILABLE = True
except Exception:
    TRIMESH_AVAILABLE = False

LOAD_CADQUERY = os.environ.get("FEF_EXPORT_STEP", "1") == "1"
try:
    if not LOAD_CADQUERY:
        raise ImportError("CadQuery import skipped because FEF_EXPORT_STEP=0")
    import cadquery as cq
    from cadquery import exporters as cq_exporters
    CADQUERY_AVAILABLE = True
except Exception:
    cq = None
    cq_exporters = None
    CADQUERY_AVAILABLE = False

# ============================================================
# 1. User configuration
# ============================================================

SEED = int(os.environ.get("FEF_SEED", "17"))
rng = np.random.default_rng(SEED)
np.random.seed(SEED)

# Moderate-precision application: desktop tracer / visual inspection / dispensing demonstrator
APPLICATION_NAME = "desktop_curve_tracing_and_low_force_dispensing_demonstrator"

# Search settings. Increase for final runs.
N_POINTS = int(os.environ.get("FEF_N_POINTS", "120"))
N_RANDOM_PER_TOPOLOGY = int(os.environ.get("FEF_N_RANDOM", "1800"))
TOPK_RANDOM_REFINE_PER_TOPOLOGY = int(os.environ.get("FEF_TOPK_REFINE", "4"))
RUN_DE_GLOBAL = os.environ.get("FEF_RUN_DE", "1") == "1"
DE_RESTARTS = int(os.environ.get("FEF_DE_RESTARTS", "1"))
DE_MAXITER = int(os.environ.get("FEF_DE_MAXITER", "45"))
DE_POPSIZE = int(os.environ.get("FEF_DE_POPSIZE", "9"))
LOCAL_MAXITER = int(os.environ.get("FEF_LOCAL_MAXITER", "260"))
SHAPE_SEEDS_PER_TOPOLOGY = int(os.environ.get("FEF_SHAPE_SEEDS", "36"))
RUN_LS_POLISH = os.environ.get("FEF_RUN_LS_POLISH", "1") == "1"
LS_MAX_NFEV = int(os.environ.get("FEF_LS_MAX_NFEV", "700"))
RANDOM_SAMPLER = os.environ.get("FEF_RANDOM_SAMPLER", "lhs").strip().lower()
ELITE_NEIGHBORHOODS = int(os.environ.get("FEF_ELITE_NEIGHBORHOODS", "8"))
ELITE_NEIGHBOR_SAMPLES = int(os.environ.get("FEF_ELITE_NEIGHBOR_SAMPLES", "6"))
ELITE_NEIGHBOR_SCALES = [
    float(v) for v in os.environ.get("FEF_ELITE_NEIGHBOR_SCALES", "0.04,0.09,0.16").split(",")
    if v.strip()
]
DIVERSE_REFINE = os.environ.get("FEF_DIVERSE_REFINE", "1") == "1"
DIVERSE_MIN_DIST = float(os.environ.get("FEF_DIVERSE_MIN_DIST", "0.08"))

# Accepted design tolerance for non-welding additive-manufacturing demonstrator.
# Use FEF_ACCEPT_RMS_MM=2.0, for example, when the current synthesis goal is
# "find a mechanism below 2 mm RMS" instead of the stricter default.
ACCEPT_RMS_MM = float(os.environ.get("FEF_ACCEPT_RMS_MM", "1.50"))
ACCEPT_MAX_MM = float(os.environ.get("FEF_ACCEPT_MAX_MM", "3.00"))
RELAXED_RMS_MM = float(os.environ.get("FEF_RELAXED_RMS_MM", str(max(3.00, 2.0 * ACCEPT_RMS_MM))))
RELAXED_MAX_MM = float(os.environ.get("FEF_RELAXED_MAX_MM", str(max(6.00, 2.0 * ACCEPT_MAX_MM))))

# Physical prototype scale
TARGET_WIDTH_MM = float(os.environ.get("FEF_TARGET_WIDTH_MM", "90.0"))
TARGET_CSV = os.environ.get("FEF_TARGET_CSV", "").strip() or None
TARGET_MODE = os.environ.get("FEF_TARGET_MODE", "hidden_stephenson").strip().lower()
GEOMETRY_PROFILE = os.environ.get("FEF_GEOMETRY_PROFILE", "compact").strip().lower()
MAX_PRINTABLE_LINK_MM = float(os.environ.get(
    "FEF_MAX_PRINTABLE_LINK_MM",
    "240.0" if GEOMETRY_PROFILE == "relaxed" else "180.0",
))
PLATE_THICKNESS_MM = 4.0
LINK_WIDTH_MM = 11.0
PIN_HOLE_DIAM_MM = 3.4      # good starting point for M3 clearance in many printers
PIN_BOSS_DIAM_MM = 7.5
BASE_THICKNESS_MM = 5.0
BASE_MARGIN_MM = 22.0

# Families to search.  The default set is intended for the Q1 capability
# benchmark: the graph/optimizer is the search engine, while these mechanism
# families are the objects being compared.
TOPOLOGIES = [
    s.strip()
    for s in os.environ.get(
        "FEF_TOPOLOGIES",
        "fourbar,slidercrank,wattI,wattII,stephensonI,stephensonII,stephensonIII,gearedFiveBar",
    ).split(",")
    if s.strip()
]

# Objective weights
W_RMS = 1.0
W_MAX = 0.22
W_SMOOTH = 0.03
W_SING = 0.08
W_SIZE = 0.015

# Output
OUT_PREFIX = os.environ.get("FEF_OUT_PREFIX", "fef_multitopology_additive_prototype")
OUT_TAG = os.environ.get("FEF_OUT_TAG", "").strip()
OUTPUT_ROOT = Path(os.environ.get("FEF_OUTPUT_ROOT", ".")).resolve()
SAVE_FIGURES = os.environ.get("FEF_SAVE_FIGURES", "1") == "1"
EXPORT_STL = os.environ.get("FEF_EXPORT_STL", "1") == "1"
EXPORT_STEP = os.environ.get("FEF_EXPORT_STEP", "1") == "1"
EXPORT_SVG = os.environ.get("FEF_EXPORT_SVG", "1") == "1"
EXPORT_ZIP = os.environ.get("FEF_EXPORT_ZIP", "1") == "1"

EPS = 1e-12
PHASE_STRIDE = int(os.environ.get("FEF_PHASE_STRIDE", "6"))
RESAMPLE_GENERATED_FOR_ERROR = os.environ.get("FEF_RESAMPLE_GENERATED_FOR_ERROR", "1") == "1"

# ============================================================
# 2. Numba kernels and geometry utilities
# ============================================================

@njit(cache=True)
def _rms_max_numba(P: np.ndarray, Y: np.ndarray) -> Tuple[float, float]:
    n = P.shape[0]
    s = 0.0
    m = 0.0
    for i in range(n):
        dx = P[i, 0] - Y[i, 0]
        dy = P[i, 1] - Y[i, 1]
        e2 = dx * dx + dy * dy
        s += e2
        e = math.sqrt(e2)
        if e > m:
            m = e
    return math.sqrt(s / n), m

@njit(cache=True)
def _smoothness_numba(P: np.ndarray) -> float:
    n = P.shape[0]
    if n < 3:
        return 0.0
    s = 0.0
    c = 0
    for i in range(1, n - 1):
        dx = P[i+1, 0] - 2.0 * P[i, 0] + P[i-1, 0]
        dy = P[i+1, 1] - 2.0 * P[i, 1] + P[i-1, 1]
        s += dx * dx + dy * dy
        c += 1
    return math.sqrt(s / max(c, 1))

def smoothness_metric(P: np.ndarray) -> float:
    if CPP_CORE_AVAILABLE and cpp_core is not None:
        return float(cpp_core.smoothness(np.ascontiguousarray(P, dtype=float)))
    return float(_smoothness_numba(P))

def rot2(theta: float) -> np.ndarray:
    c, s = math.cos(theta), math.sin(theta)
    return np.array([[c, -s], [s, c]], dtype=float)

def norm2(v: np.ndarray) -> float:
    return float(np.linalg.norm(v))

def unit(v: np.ndarray) -> np.ndarray:
    n = norm2(v)
    if n < EPS:
        raise ValueError("zero vector")
    return v / n

def perp(v: np.ndarray) -> np.ndarray:
    return np.array([-v[1], v[0]], dtype=float)

def cross2(a: np.ndarray, b: np.ndarray) -> float:
    return float(a[0] * b[1] - a[1] * b[0])

def circle_intersections(c0: np.ndarray, r0: float, c1: np.ndarray, r1: float) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    dvec = c1 - c0
    d = norm2(dvec)
    if d < 1e-12:
        return None
    if d > r0 + r1 + 1e-9:
        return None
    if d < abs(r0 - r1) - 1e-9:
        return None
    a = (r0*r0 - r1*r1 + d*d) / (2.0*d)
    h2 = r0*r0 - a*a
    if h2 < -1e-9:
        return None
    h = math.sqrt(max(h2, 0.0))
    ex = dvec / d
    p2 = c0 + a * ex
    ey = np.array([-ex[1], ex[0]], dtype=float)
    return p2 + h * ey, p2 - h * ey

def choose_branch(cands: Tuple[np.ndarray, np.ndarray], prev: Optional[np.ndarray], prefer_positive_y: bool = True) -> np.ndarray:
    p1, p2 = cands
    if prev is None:
        if prefer_positive_y:
            return p1 if p1[1] >= p2[1] else p2
        return p1
    return p1 if norm2(p1 - prev) <= norm2(p2 - prev) else p2

def procrustes_rigid(X: np.ndarray, Y: np.ndarray) -> Tuple[np.ndarray, float, np.ndarray, np.ndarray]:
    Xc = X - X.mean(axis=0)
    Yc = Y - Y.mean(axis=0)
    H = Xc.T @ Yc
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1.0
        R = Vt.T @ U.T
    t = Y.mean(axis=0) - (R @ X.mean(axis=0))
    Xa = (R @ X.T).T + t
    rms, _ = _rms_max_numba(Xa, Y)
    return Xa, float(rms), R, t

def best_alignment(P: np.ndarray, Y: np.ndarray, allow_phase: bool = True, phase_stride: int = PHASE_STRIDE) -> Dict[str, Any]:
    """Rigid alignment + cyclic phase. Uses a coarse-to-local phase search for speed."""
    if CPP_CORE_AVAILABLE and cpp_core is not None:
        return cpp_core.best_alignment(
            np.ascontiguousarray(P, dtype=float),
            np.ascontiguousarray(Y, dtype=float),
            bool(allow_phase),
            int(phase_stride),
        )
    n = P.shape[0]
    if not allow_phase:
        candidate_shifts = [0]
    else:
        stride = max(1, int(phase_stride))
        coarse = list(range(0, n, stride))
        best_s = 0
        best_r = np.inf
        # Coarse phase scan
        for s in coarse:
            Ys = np.roll(Y, shift=s, axis=0)
            Pa, rms, R, t = procrustes_rigid(P, Ys)
            if rms < best_r:
                best_r = rms
                best_s = s
        # Local refinement around coarse best
        rad = max(1, stride)
        candidate_shifts = sorted(set([(best_s + ds) % n for ds in range(-rad, rad + 1)]))

    best = {"rms": np.inf, "max": np.inf, "shift": 0, "aligned": None, "target_shifted": None}
    for s in candidate_shifts:
        Ys = np.roll(Y, shift=s, axis=0)
        Pa, rms, R, t = procrustes_rigid(P, Ys)
        rms2, maxe = _rms_max_numba(Pa, Ys)
        if rms2 < best["rms"]:
                best = {"rms": float(rms2), "max": float(maxe), "shift": int(s), "aligned": Pa, "target_shifted": Ys, "R": R, "t": t}
    return best

def alignment_at_shift(P: np.ndarray, Y: np.ndarray, shift: int) -> Dict[str, Any]:
    if CPP_CORE_AVAILABLE and cpp_core is not None:
        return cpp_core.alignment_at_shift(
            np.ascontiguousarray(P, dtype=float),
            np.ascontiguousarray(Y, dtype=float),
            int(shift),
        )
    Ys = np.roll(Y, shift=int(shift), axis=0)
    Pa, rms, R, t = procrustes_rigid(P, Ys)
    rms2, maxe = _rms_max_numba(Pa, Ys)
    return {
        "rms": float(rms2),
        "max": float(maxe),
        "shift": int(shift),
        "aligned": Pa,
        "target_shifted": Ys,
        "R": R,
        "t": t,
    }

def path_width(P: np.ndarray) -> float:
    return float(np.max(P[:, 0]) - np.min(P[:, 0]))

def path_height(P: np.ndarray) -> float:
    return float(np.max(P[:, 1]) - np.min(P[:, 1]))

def normalize_target_width(P: np.ndarray, width_mm: float) -> np.ndarray:
    P0 = P - P.mean(axis=0)
    w = path_width(P0)
    if w < EPS:
        return P0
    return P0 * (width_mm / w)

def resample_closed_curve(P: np.ndarray, n: int) -> np.ndarray:
    """Resample an ordered closed curve to n approximately arclength-uniform points."""
    if CPP_CORE_AVAILABLE and cpp_core is not None:
        return cpp_core.resample_closed_curve(np.ascontiguousarray(P, dtype=float), int(n))
    P = np.asarray(P, dtype=float)
    if P.ndim != 2 or P.shape[1] != 2 or len(P) < 4:
        raise ValueError("target curve must be an Nx2 array with at least four points")
    if norm2(P[0] - P[-1]) < 1e-9:
        P = P[:-1]
    closed = np.vstack([P, P[0]])
    seg = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    s = np.concatenate([[0.0], np.cumsum(seg)])
    total = float(s[-1])
    if total < EPS:
        raise ValueError("target curve has near-zero arclength")
    sample_s = np.linspace(0.0, total, n + 1, endpoint=True)[:-1]
    x = np.interp(sample_s, s, closed[:, 0])
    y = np.interp(sample_s, s, closed[:, 1])
    return np.column_stack([x, y])

def curve_for_error(P: np.ndarray, n: int) -> np.ndarray:
    """Return the mechanism path used by the geometric error metric.

    Mechanism points are naturally sampled by input angle, while a target curve is
    usually a geometric path without timing.  Resampling the generated path by
    arclength keeps the benchmark focused on shape synthesis instead of input
    speed matching.
    """
    if not RESAMPLE_GENERATED_FOR_ERROR:
        return P
    try:
        return resample_closed_curve(P, n)
    except Exception:
        return P

def curve_complexity_metrics(P: np.ndarray) -> Dict[str, float]:
    P = np.asarray(P, dtype=float)
    closed = np.vstack([P, P[0]])
    seg = np.linalg.norm(np.diff(closed, axis=0), axis=1)
    perimeter = float(np.sum(seg))
    chord = float(max(path_width(P), path_height(P), EPS))
    d1 = np.roll(P, -1, axis=0) - np.roll(P, 1, axis=0)
    d2 = np.roll(P, -1, axis=0) - 2.0 * P + np.roll(P, 1, axis=0)
    cross = d1[:, 0] * d2[:, 1] - d1[:, 1] * d2[:, 0]
    speed = np.linalg.norm(d1, axis=1)
    curvature = np.abs(cross) / np.maximum(speed ** 3, EPS)
    turn = np.arctan2(cross, np.sum(d1 * np.roll(d1, -1, axis=0), axis=1))
    sign_changes = int(np.sum(np.diff(np.signbit(cross)).astype(int) != 0))
    return {
        "perimeter_mm": perimeter,
        "perimeter_to_width": float(perimeter / chord),
        "curvature_mean": float(np.mean(curvature)),
        "curvature_p95": float(np.percentile(curvature, 95)),
        "turning_energy": float(np.mean(turn ** 2)),
        "curvature_sign_changes": float(sign_changes),
    }

def generate_independent_benchmark_target(theta: np.ndarray, difficulty: str) -> Tuple[np.ndarray, Dict[str, Any]]:
    t = np.asarray(theta, dtype=float)
    difficulty = difficulty.lower()
    if difficulty == "easy":
        x = np.cos(t) + 0.10 * np.cos(2.0 * t + 0.25)
        y = 0.55 * np.sin(t) + 0.05 * np.sin(2.0 * t - 0.40)
        description = "smooth oval with weak second harmonic"
    elif difficulty == "medium":
        x = np.cos(t) + 0.22 * np.cos(2.0 * t - 0.35) - 0.08 * np.cos(3.0 * t + 0.60)
        y = 0.68 * np.sin(t) + 0.20 * np.sin(2.0 * t + 0.25) + 0.08 * np.sin(3.0 * t - 0.80)
        description = "asymmetric closed curve with moderate harmonic content"
    elif difficulty in {"medium_hard", "moderate_hard"}:
        xm = np.cos(t) + 0.22 * np.cos(2.0 * t - 0.35) - 0.08 * np.cos(3.0 * t + 0.60)
        ym = 0.68 * np.sin(t) + 0.20 * np.sin(2.0 * t + 0.25) + 0.08 * np.sin(3.0 * t - 0.80)
        xh = np.cos(t) + 0.36 * np.cos(2.0 * t + 0.20) - 0.24 * np.cos(3.0 * t - 0.65) + 0.11 * np.cos(5.0 * t)
        yh = 0.74 * np.sin(t) - 0.34 * np.sin(2.0 * t - 0.30) + 0.20 * np.sin(4.0 * t + 0.45) - 0.08 * np.sin(5.0 * t)
        alpha = 0.35
        x = (1.0 - alpha) * xm + alpha * xh
        y = (1.0 - alpha) * ym + alpha * yh
        difficulty = "medium_hard"
        description = "35 percent interpolation between the moderate and high-curvature benchmarks"
    elif difficulty == "hard":
        x = np.cos(t) + 0.36 * np.cos(2.0 * t + 0.20) - 0.24 * np.cos(3.0 * t - 0.65) + 0.11 * np.cos(5.0 * t)
        y = 0.74 * np.sin(t) - 0.34 * np.sin(2.0 * t - 0.30) + 0.20 * np.sin(4.0 * t + 0.45) - 0.08 * np.sin(5.0 * t)
        description = "high-curvature multi-harmonic curve"
    else:
        raise ValueError(f"Unknown FEF_TARGET_MODE={difficulty!r}; use easy, medium, medium_hard, hard, hidden_stephenson, or FEF_TARGET_CSV.")

    raw = np.column_stack([x, y]).astype(float)
    target = normalize_target_width(resample_closed_curve(raw, len(theta)), TARGET_WIDTH_MM)
    metrics = curve_complexity_metrics(target)
    return target, {
        "target_source": f"independent_fourier_{difficulty}",
        "target_mode": difficulty,
        "generated_by_stephenson": False,
        "application": APPLICATION_NAME,
        "target_width_mm": TARGET_WIDTH_MM,
        "description": description,
        "n_resampled_points": int(len(theta)),
        "complexity_metrics": metrics,
        "tolerance_rms_mm": ACCEPT_RMS_MM,
        "tolerance_max_mm": ACCEPT_MAX_MM,
    }

def load_or_generate_target(theta: np.ndarray) -> Tuple[np.ndarray, Dict[str, Any]]:
    if TARGET_CSV:
        df = pd.read_csv(TARGET_CSV)
        if not {"x", "y"}.issubset(df.columns):
            raise ValueError("FEF_TARGET_CSV must contain x and y columns")
        raw = df[["x", "y"]].to_numpy(dtype=float)
        target = resample_closed_curve(raw, len(theta))
        target = normalize_target_width(target, TARGET_WIDTH_MM)
        return target, {
            "target_source": "csv_curve",
            "target_csv": str(TARGET_CSV),
            "generated_by_stephenson": False,
            "application": APPLICATION_NAME,
            "target_width_mm": TARGET_WIDTH_MM,
            "n_input_points": int(len(raw)),
            "n_resampled_points": int(len(theta)),
            "tolerance_rms_mm": ACCEPT_RMS_MM,
            "tolerance_max_mm": ACCEPT_MAX_MM,
        }
    if TARGET_MODE in {"easy", "medium", "medium_hard", "moderate_hard", "hard"}:
        return generate_independent_benchmark_target(theta, TARGET_MODE)
    if TARGET_MODE not in {"hidden_stephenson", "application", "default"}:
        raise ValueError(f"Unknown FEF_TARGET_MODE={TARGET_MODE!r}; use easy, medium, medium_hard, hard, hidden_stephenson, or FEF_TARGET_CSV.")
    return generate_application_target(theta)

# ============================================================
# 3. Targets and mechanism models
# ============================================================

def theta_samples(n: int = N_POINTS) -> np.ndarray:
    return np.linspace(0.0, 2.0 * math.pi, n, endpoint=False)

@dataclass
class SimResult:
    success: bool
    P: Optional[np.ndarray]
    data: Optional[Dict[str, np.ndarray]]
    min_delta: float = 0.0
    message: str = ""

# ---------- Four-bar ----------
@dataclass
class FourBarParams:
    L1: float
    L2: float
    L3: float
    L4: float
    rP: float
    phiP: float

def simulate_fourbar_params(p: FourBarParams, theta: np.ndarray) -> SimResult:
    O2 = np.array([0.0, 0.0])
    O4 = np.array([p.L1, 0.0])
    A_all, B_all, P_all = [], [], []
    prev_B = None
    deltas = []
    for th in theta:
        A = O2 + p.L2 * np.array([math.cos(th), math.sin(th)])
        cand = circle_intersections(A, p.L3, O4, p.L4)
        if cand is None:
            return SimResult(False, None, None, message="closure")
        B = choose_branch(cand, prev_B, prefer_positive_y=True)
        u = unit(B - A)
        v = perp(u)
        P = A + p.rP * (math.cos(p.phiP) * u + math.sin(p.phiP) * v)
        delta = abs(cross2(B - A, O4 - B))
        A_all.append(A); B_all.append(B); P_all.append(P); deltas.append(delta)
        prev_B = B.copy()
    data = {"O2": O2, "O4": O4, "A": np.asarray(A_all), "B": np.asarray(B_all)}
    return SimResult(True, np.asarray(P_all), data, float(np.min(deltas)))

# ---------- Slider-crank ----------
@dataclass
class SliderCrankParams:
    r: float
    l: float
    q: float
    h: float
    x0: float

def simulate_slider_params(p: SliderCrankParams, theta: np.ndarray) -> SimResult:
    if p.l <= p.r:
        return SimResult(False, None, None, message="l<=r")
    O = np.array([0.0, 0.0])
    A_all, B_all, P_all = [], [], []
    for th in theta:
        A = O + p.r * np.array([math.cos(th), math.sin(th)])
        rad = p.l * p.l - A[1] * A[1]
        if rad < 0.0:
            return SimResult(False, None, None, message="slider closure")
        B = np.array([p.x0 + A[0] + math.sqrt(max(rad, 0.0)), 0.0])
        u = unit(B - A)
        v = perp(u)
        P = A + p.q * p.l * u + p.h * v
        A_all.append(A); B_all.append(B); P_all.append(P)
    data = {"O": O, "A": np.asarray(A_all), "B": np.asarray(B_all)}
    # crude margin: connecting rod vertical feasibility
    return SimResult(True, np.asarray(P_all), data, float(p.l - p.r))

# ---------- Stephenson III exact ----------
@dataclass
class StephensonIIIParams:
    L0: float
    L2: float
    L3: float
    L4: float
    r_c: float
    phi_c: float
    L5: float
    L6: float
    pu: float
    pv: float

def simulate_stephenson_params(p: StephensonIIIParams, theta: np.ndarray) -> SimResult:
    O2 = np.array([0.0, 0.0])
    O4 = np.array([p.L0, 0.0])
    A_all, B_all, C_all, D_all, P_all = [], [], [], [], []
    prev_B = None
    prev_D = None
    deltas = []
    for th in theta:
        A = O2 + p.L2 * np.array([math.cos(th), math.sin(th)])
        cand_B = circle_intersections(A, p.L3, O4, p.L4)
        if cand_B is None:
            return SimResult(False, None, None, message="loop1")
        B = choose_branch(cand_B, prev_B, prefer_positive_y=True)
        uAB = unit(B - A)
        vAB = perp(uAB)
        C = A + p.r_c * (math.cos(p.phi_c) * uAB + math.sin(p.phi_c) * vAB)
        cand_D = circle_intersections(C, p.L5, O4, p.L6)
        if cand_D is None:
            return SimResult(False, None, None, message="loop2")
        D = choose_branch(cand_D, prev_D, prefer_positive_y=False)
        uCD = unit(D - C)
        vCD = perp(uCD)
        # Output point on CD link: this is the useful Stephenson output for the new prototype study.
        P = C + p.pu * uCD + p.pv * vCD
        delta1 = abs(cross2(B - A, O4 - B))
        delta2 = abs(cross2(D - C, O4 - D))
        deltas.append(min(delta1, delta2))
        A_all.append(A); B_all.append(B); C_all.append(C); D_all.append(D); P_all.append(P)
        prev_B = B.copy(); prev_D = D.copy()
    data = {"O2": O2, "O4": O4, "A": np.asarray(A_all), "B": np.asarray(B_all), "C": np.asarray(C_all), "D": np.asarray(D_all)}
    return SimResult(True, np.asarray(P_all), data, float(np.min(deltas)))

# ---------- Watt and Stephenson variants for family-capability benchmarks ----------
@dataclass
class WattIParams:
    L0: float
    L2: float
    L3: float
    L4: float
    r_c: float
    phi_c: float
    gx: float
    gy: float
    L5: float
    L6: float
    pu: float
    pv: float

@dataclass
class StephensonDyadParams:
    L0: float
    L2: float
    L3: float
    L4: float
    L5: float
    L6: float
    pu: float
    pv: float

@dataclass
class GearedFiveBarParams:
    L0: float
    L2: float
    L4: float
    L3: float
    L5: float
    ratio: float
    phase: float
    r_p: float
    phi_p: float

def _fourbar_loop_points(L0: float, L2: float, L3: float, L4: float, theta: np.ndarray):
    O2 = np.array([0.0, 0.0])
    O4 = np.array([L0, 0.0])
    A_all, B_all, deltas = [], [], []
    prev_B = None
    for th in theta:
        A = O2 + L2 * np.array([math.cos(th), math.sin(th)])
        cand_B = circle_intersections(A, L3, O4, L4)
        if cand_B is None:
            return False, O2, O4, None, None, 0.0
        B = choose_branch(cand_B, prev_B, prefer_positive_y=True)
        A_all.append(A)
        B_all.append(B)
        deltas.append(abs(cross2(B - A, O4 - B)))
        prev_B = B.copy()
    return True, O2, O4, np.asarray(A_all), np.asarray(B_all), float(np.min(deltas))

def _simulate_grounded_second_dyad(
    theta: np.ndarray,
    O2: np.ndarray,
    O4: np.ndarray,
    A_all: np.ndarray,
    B_all: np.ndarray,
    anchors: np.ndarray,
    ground: np.ndarray,
    L5: float,
    L6: float,
    pu: float,
    pv: float,
    prefer_positive_y: bool = False,
) -> SimResult:
    D_all, P_all, deltas = [], [], []
    prev_D = None
    for C in anchors:
        cand_D = circle_intersections(C, L5, ground, L6)
        if cand_D is None:
            return SimResult(False, None, None, message="second-loop closure")
        D = choose_branch(cand_D, prev_D, prefer_positive_y=prefer_positive_y)
        u = unit(D - C)
        v = perp(u)
        P = C + pu * u + pv * v
        deltas.append(abs(cross2(D - C, ground - D)))
        D_all.append(D)
        P_all.append(P)
        prev_D = D.copy()
    data = {"O2": O2, "O4": O4, "A": A_all, "B": B_all, "C": anchors, "D": np.asarray(D_all)}
    return SimResult(True, np.asarray(P_all), data, float(np.min(deltas)))

def simulate_wattI_params(p: WattIParams, theta: np.ndarray) -> SimResult:
    ok, O2, O4, A_all, B_all, min_delta1 = _fourbar_loop_points(p.L0, p.L2, p.L3, p.L4, theta)
    if not ok:
        return SimResult(False, None, None, message="loop1")
    C_all = []
    for A, B in zip(A_all, B_all):
        uAB = unit(B - A)
        vAB = perp(uAB)
        C_all.append(A + p.r_c * (math.cos(p.phi_c) * uAB + math.sin(p.phi_c) * vAB))
    sim = _simulate_grounded_second_dyad(theta, O2, O4, A_all, B_all, np.asarray(C_all), np.array([p.gx, p.gy]), p.L5, p.L6, p.pu, p.pv)
    sim.min_delta = min(float(sim.min_delta), min_delta1) if sim.success else 0.0
    return sim

def simulate_wattII_params(p: WattIParams, theta: np.ndarray) -> SimResult:
    ok, O2, O4, A_all, B_all, min_delta1 = _fourbar_loop_points(p.L0, p.L2, p.L3, p.L4, theta)
    if not ok:
        return SimResult(False, None, None, message="loop1")
    C_all = []
    for B in B_all:
        uBO4 = unit(O4 - B)
        vBO4 = perp(uBO4)
        C_all.append(B + p.r_c * (math.cos(p.phi_c) * uBO4 + math.sin(p.phi_c) * vBO4))
    sim = _simulate_grounded_second_dyad(theta, O2, O4, A_all, B_all, np.asarray(C_all), np.array([p.gx, p.gy]), p.L5, p.L6, p.pu, p.pv)
    sim.min_delta = min(float(sim.min_delta), min_delta1) if sim.success else 0.0
    return sim

def simulate_stephensonI_params(p: StephensonDyadParams, theta: np.ndarray) -> SimResult:
    ok, O2, O4, A_all, B_all, min_delta1 = _fourbar_loop_points(p.L0, p.L2, p.L3, p.L4, theta)
    if not ok:
        return SimResult(False, None, None, message="loop1")
    sim = _simulate_grounded_second_dyad(theta, O2, O4, A_all, B_all, A_all, O4, p.L5, p.L6, p.pu, p.pv)
    sim.min_delta = min(float(sim.min_delta), min_delta1) if sim.success else 0.0
    return sim

def simulate_stephensonII_params(p: StephensonDyadParams, theta: np.ndarray) -> SimResult:
    ok, O2, O4, A_all, B_all, min_delta1 = _fourbar_loop_points(p.L0, p.L2, p.L3, p.L4, theta)
    if not ok:
        return SimResult(False, None, None, message="loop1")
    sim = _simulate_grounded_second_dyad(theta, O2, O4, A_all, B_all, B_all, O2, p.L5, p.L6, p.pu, p.pv, prefer_positive_y=True)
    sim.min_delta = min(float(sim.min_delta), min_delta1) if sim.success else 0.0
    return sim

def simulate_geared_fivebar_params(p: GearedFiveBarParams, theta: np.ndarray) -> SimResult:
    O2 = np.array([0.0, 0.0])
    O4 = np.array([p.L0, 0.0])
    A_all, B_all, C_all, P_all, deltas = [], [], [], [], []
    prev_C = None
    for th in theta:
        A = O2 + p.L2 * np.array([math.cos(th), math.sin(th)])
        B = O4 + p.L4 * np.array([math.cos(p.ratio * th + p.phase), math.sin(p.ratio * th + p.phase)])
        cand_C = circle_intersections(A, p.L3, B, p.L5)
        if cand_C is None:
            return SimResult(False, None, None, message="fivebar closure")
        C = choose_branch(cand_C, prev_C, prefer_positive_y=True)
        u = unit(C - A)
        v = perp(u)
        P = C + p.r_p * (math.cos(p.phi_p) * u + math.sin(p.phi_p) * v)
        deltas.append(min(abs(cross2(C - A, B - C)), abs(cross2(A - O2, B - O4))))
        A_all.append(A); B_all.append(B); C_all.append(C); P_all.append(P)
        prev_C = C.copy()
    data = {"O2": O2, "O4": O4, "A": np.asarray(A_all), "B": np.asarray(B_all), "C": np.asarray(C_all)}
    return SimResult(True, np.asarray(P_all), data, float(np.min(deltas)))

# ---------- parameter definitions ----------
PARAMS: Dict[str, Dict[str, Any]] = {
    "fourbar": {
        "names": ["L1", "L2", "L3", "L4", "rP", "phiP"],
        "lb": np.array([35.0, 12.0, 20.0, 20.0, 0.0, -math.pi]),
        "ub": np.array([150.0, 90.0, 165.0, 165.0, 150.0, math.pi]),
    },
    "slidercrank": {
        "names": ["r", "l", "q", "h", "x0"],
        "lb": np.array([8.0, 25.0, -0.30, -60.0, -30.0]),
        "ub": np.array([80.0, 210.0, 1.35, 60.0, 60.0]),
    },
    "wattI": {
        "names": ["L0", "L2", "L3", "L4", "r_c", "phi_c", "gx", "gy", "L5", "L6", "pu", "pv"],
        "lb": np.array([35.0, 12.0, 20.0, 20.0, 0.0, -math.pi, -40.0, -90.0, 20.0, 20.0, -50.0, -80.0]),
        "ub": np.array([160.0, 95.0, 180.0, 180.0, 150.0, math.pi, 210.0, 90.0, 180.0, 180.0, 180.0, 80.0]),
    },
    "wattII": {
        "names": ["L0", "L2", "L3", "L4", "r_c", "phi_c", "gx", "gy", "L5", "L6", "pu", "pv"],
        "lb": np.array([35.0, 12.0, 20.0, 20.0, 0.0, -math.pi, -40.0, -90.0, 20.0, 20.0, -50.0, -80.0]),
        "ub": np.array([160.0, 95.0, 180.0, 180.0, 150.0, math.pi, 210.0, 90.0, 180.0, 180.0, 180.0, 80.0]),
    },
    "stephensonI": {
        "names": ["L0", "L2", "L3", "L4", "L5", "L6", "pu", "pv"],
        "lb": np.array([35.0, 12.0, 20.0, 20.0, 20.0, 20.0, -50.0, -80.0]),
        "ub": np.array([160.0, 95.0, 180.0, 180.0, 180.0, 180.0, 180.0, 80.0]),
    },
    "stephensonII": {
        "names": ["L0", "L2", "L3", "L4", "L5", "L6", "pu", "pv"],
        "lb": np.array([35.0, 12.0, 20.0, 20.0, 20.0, 20.0, -50.0, -80.0]),
        "ub": np.array([160.0, 95.0, 180.0, 180.0, 180.0, 180.0, 180.0, 80.0]),
    },
    "stephensonIII": {
        "names": ["L0", "L2", "L3", "L4", "r_c", "phi_c", "L5", "L6", "pu", "pv"],
        "lb": np.array([35.0, 12.0, 20.0, 20.0, 0.0, -math.pi, 20.0, 20.0, -30.0, -70.0]),
        "ub": np.array([160.0, 95.0, 180.0, 180.0, 150.0, math.pi, 180.0, 180.0, 180.0, 70.0]),
    },
    "gearedFiveBar": {
        "names": ["L0", "L2", "L4", "L3", "L5", "ratio", "phase", "r_p", "phi_p"],
        "lb": np.array([35.0, 10.0, 10.0, 20.0, 20.0, -3.0, -math.pi, 0.0, -math.pi]),
        "ub": np.array([160.0, 95.0, 95.0, 180.0, 180.0, 3.0, math.pi, 150.0, math.pi]),
    },
}

def apply_geometry_profile(params: Dict[str, Dict[str, Any]]) -> None:
    """Adjust search bounds for compact or relaxed publication experiments.

    The compact profile is intended for a desktop printable demonstrator.  The
    relaxed profile keeps the same mechanism families but allows a larger
    laboratory-scale prototype, which is useful for testing whether failure is
    caused by kinematic family limits or by overly tight packaging constraints.
    """
    if GEOMETRY_PROFILE in {"compact", ""}:
        return
    if GEOMETRY_PROFILE != "relaxed":
        raise ValueError("FEF_GEOMETRY_PROFILE must be compact or relaxed")

    link_max = max(240.0, MAX_PRINTABLE_LINK_MM)
    ground_max = max(220.0, 0.90 * link_max)
    crank_max = max(140.0, 0.60 * link_max)
    marker_max = max(260.0, 1.10 * link_max)
    vertical_max = max(140.0, 0.60 * link_max)
    external_ground_x = max(300.0, 1.25 * link_max)
    external_ground_y = max(160.0, 0.70 * link_max)
    relaxed_ub = {
        "L0": ground_max, "L1": ground_max,
        "L2": crank_max, "r": crank_max,
        "L3": link_max, "L4": link_max, "L5": link_max, "L6": link_max,
        "l": 300.0,
        "rP": marker_max, "r_c": marker_max, "r_p": marker_max,
        "gx": external_ground_x, "gy": external_ground_y,
        "pu": marker_max, "pv": vertical_max,
        "h": 100.0, "x0": 100.0,
    }
    relaxed_lb = {
        "gx": -100.0, "gy": -external_ground_y,
        "pu": -0.50 * marker_max, "pv": -vertical_max,
        "h": -100.0, "x0": -70.0,
    }
    for spec in params.values():
        names = spec["names"]
        lb = spec["lb"].copy()
        ub = spec["ub"].copy()
        for i, name in enumerate(names):
            if name in relaxed_ub:
                ub[i] = max(float(ub[i]), relaxed_ub[name])
            if name in relaxed_lb:
                lb[i] = min(float(lb[i]), relaxed_lb[name])
        spec["lb"] = lb
        spec["ub"] = ub

apply_geometry_profile(PARAMS)

def vector_to_params(topo: str, x: np.ndarray):
    x = np.asarray(x, dtype=float)
    names = PARAMS[topo]["names"]
    d = {k: float(v) for k, v in zip(names, x)}
    if topo == "fourbar":
        return FourBarParams(**d)
    if topo == "slidercrank":
        return SliderCrankParams(**d)
    if topo in {"wattI", "wattII"}:
        return WattIParams(**d)
    if topo in {"stephensonI", "stephensonII"}:
        return StephensonDyadParams(**d)
    if topo == "stephensonIII":
        return StephensonIIIParams(**d)
    if topo == "gearedFiveBar":
        return GearedFiveBarParams(**d)
    raise ValueError(topo)

def simulate_topology(topo: str, x: np.ndarray, theta: np.ndarray) -> SimResult:
    try:
        p = vector_to_params(topo, x)
        if topo == "fourbar":
            return simulate_fourbar_params(p, theta)
        if topo == "slidercrank":
            return simulate_slider_params(p, theta)
        if topo == "wattI":
            return simulate_wattI_params(p, theta)
        if topo == "wattII":
            return simulate_wattII_params(p, theta)
        if topo == "stephensonI":
            return simulate_stephensonI_params(p, theta)
        if topo == "stephensonII":
            return simulate_stephensonII_params(p, theta)
        if topo == "stephensonIII":
            return simulate_stephenson_params(p, theta)
        if topo == "gearedFiveBar":
            return simulate_geared_fivebar_params(p, theta)
    except Exception as e:
        return SimResult(False, None, None, message=str(e))
    return SimResult(False, None, None, message="unknown")

# ============================================================
# 4. Application target
# ============================================================

def generate_application_target(theta: np.ndarray) -> Tuple[np.ndarray, Dict[str, Any]]:
    """A compatible but non-trivial target for a low-force desktop tracer.

    It is generated by a hidden Stephenson III, transformed and normalized to a physical width.
    The hidden source is saved only as metadata for audit; the optimizer does not use it.
    """
    true = StephensonIIIParams(
        L0=92.0,
        L2=34.0,
        L3=78.0,
        L4=63.0,
        r_c=46.0,
        phi_c=0.88,
        L5=52.0,
        L6=79.0,
        pu=34.0,
        pv=-18.0,
    )
    sim = simulate_stephenson_params(true, theta)
    if not sim.success:
        raise RuntimeError("Internal target Stephenson parameters are infeasible.")
    P = sim.P.copy()
    R = rot2(math.radians(19.0))
    P = (R @ P.T).T + np.array([33.0, -18.0])
    P = np.roll(P, 17, axis=0)
    P = normalize_target_width(P, TARGET_WIDTH_MM)
    meta = {
        "target_source": "hidden_stephensonIII_for_desktop_curve_tracer",
        "generated_by_stephenson": True,
        "application": APPLICATION_NAME,
        "target_width_mm": TARGET_WIDTH_MM,
        "hidden_params": asdict(true),
        "tolerance_rms_mm": ACCEPT_RMS_MM,
        "tolerance_max_mm": ACCEPT_MAX_MM,
    }
    return P, meta

# ============================================================
# 5. Evaluation and optimization
# ============================================================

def evaluate_candidate(topo: str, x: np.ndarray, target: np.ndarray, theta: np.ndarray, allow_phase: bool = True) -> Dict[str, Any]:
    lb, ub = PARAMS[topo]["lb"], PARAMS[topo]["ub"]
    x = np.minimum(np.maximum(np.asarray(x, dtype=float), lb), ub)
    sim = simulate_topology(topo, x, theta)
    if not sim.success or sim.P is None:
        return {"valid": False, "cost": 1e9, "message": sim.message, "topology": topo, "x": x}
    P_eval = curve_for_error(sim.P, len(target))
    align = best_alignment(P_eval, target, allow_phase=allow_phase)
    sm = smoothness_metric(align["aligned"])
    size = float(np.sum(x[: min(6, len(x))]))
    sing_pen = 0.0 if sim.min_delta > 50.0 else (50.0 - sim.min_delta) / 50.0
    size_pen = max(0.0, size - 650.0) / 650.0
    cost = W_RMS * align["rms"] + W_MAX * align["max"] + W_SMOOTH * sm + W_SING * sing_pen + W_SIZE * size_pen
    accepted = bool(align["rms"] <= ACCEPT_RMS_MM and align["max"] <= ACCEPT_MAX_MM)
    relaxed = bool(align["rms"] <= RELAXED_RMS_MM and align["max"] <= RELAXED_MAX_MM)
    return {
        "valid": True,
        "accepted": accepted,
        "relaxed": relaxed,
        "cost": float(cost),
        "rms_mm": float(align["rms"]),
        "max_mm": float(align["max"]),
        "smoothness": float(sm),
        "min_delta": float(sim.min_delta),
        "shift": int(align["shift"]),
        "topology": topo,
        "x": x,
        "P_raw": sim.P,
        "P_eval": P_eval,
        "P_aligned": align["aligned"],
        "target_shifted": align["target_shifted"],
        "data": sim.data,
    }

def objective_for_topology(topo: str, target: np.ndarray, theta: np.ndarray):
    def f(z):
        return evaluate_candidate(topo, z, target, theta, allow_phase=True)["cost"]
    return f

def residual_for_topology_ls(topo: str, x: np.ndarray, target: np.ndarray, theta: np.ndarray, fixed_shift: int) -> np.ndarray:
    lb, ub = PARAMS[topo]["lb"], PARAMS[topo]["ub"]
    x = np.minimum(np.maximum(np.asarray(x, dtype=float), lb), ub)
    sim = simulate_topology(topo, x, theta)
    if not sim.success or sim.P is None:
        return np.ones(2 * len(theta) + 3, dtype=float) * 1e3
    P_eval = curve_for_error(sim.P, len(target))
    align = alignment_at_shift(P_eval, target, fixed_shift)
    res = (align["aligned"] - align["target_shifted"]).reshape(-1)
    sing_pen = max(0.0, 50.0 - float(sim.min_delta)) / 50.0
    max_pen = max(0.0, float(align["max"]) - ACCEPT_MAX_MM) / max(ACCEPT_MAX_MM, EPS)
    size = float(np.sum(x[: min(6, len(x))]))
    size_pen = max(0.0, size - 650.0) / 650.0
    return np.concatenate([res, np.array([2.0 * max_pen, sing_pen, 0.25 * size_pen], dtype=float)])

def polish_candidate_with_ls(topo: str, sol: Dict[str, Any], target: np.ndarray, theta: np.ndarray) -> Dict[str, Any]:
    if not RUN_LS_POLISH or not sol.get("valid", False):
        return sol
    lb, ub = PARAMS[topo]["lb"], PARAMS[topo]["ub"]
    best = sol
    messages: List[str] = []
    nfev = 0
    shift = int(best.get("shift", 0))
    for round_id in range(2):
        try:
            ls = least_squares(
                lambda z: residual_for_topology_ls(topo, z, target, theta, fixed_shift=shift),
                np.minimum(np.maximum(best["x"], lb), ub),
                bounds=(lb, ub),
                method="trf",
                loss="soft_l1",
                f_scale=1.0,
                x_scale=np.maximum(ub - lb, 1.0),
                max_nfev=max(25, LS_MAX_NFEV // 2),
                verbose=0,
            )
            nfev += int(getattr(ls, "nfev", 0))
            ev_ls = evaluate_candidate(topo, ls.x, target, theta, allow_phase=True)
            messages.append(f"ls{round_id + 1}_rms={ev_ls.get('rms_mm', np.inf):.4f}")
            if ev_ls.get("valid", False) and ev_ls["cost"] < best["cost"]:
                best = ev_ls
                shift = int(best.get("shift", shift))
            else:
                break
        except Exception as exc:
            messages.append(f"ls_exception={type(exc).__name__}:{exc}")
            break
    best["ls_polish_nfev"] = int(nfev)
    best["ls_polish_trace"] = " | ".join(messages)
    return best

def refine_candidate(topo: str, x0: np.ndarray, target: np.ndarray, theta: np.ndarray) -> Dict[str, Any]:
    f = objective_for_topology(topo, target, theta)
    lb, ub = PARAMS[topo]["lb"], PARAMS[topo]["ub"]
    bounds = list(zip(lb, ub))
    before = evaluate_candidate(topo, x0, target, theta)
    # Powell is robust for the nonsmooth phase/alignment objective.
    res = minimize(
        f,
        np.minimum(np.maximum(x0, lb), ub),
        method="Powell",
        bounds=bounds,
        options={"maxiter": LOCAL_MAXITER, "xtol": 1e-4, "ftol": 1e-4, "disp": False},
    )
    after = evaluate_candidate(topo, res.x, target, theta)
    # Keep the better one to avoid degradation.
    best = after if after["cost"] < before["cost"] else before
    best = polish_candidate_with_ls(topo, best, target, theta)
    best["refine_success"] = bool(res.success)
    best["refine_nfev"] = int(getattr(res, "nfev", 0)) + int(best.get("ls_polish_nfev", 0))
    best["refine_message"] = str(res.message)
    best["before_rms_mm"] = float(before.get("rms_mm", np.nan))
    best["before_max_mm"] = float(before.get("max_mm", np.nan))
    return best

def random_sample_topology(topo: str, n: int) -> np.ndarray:
    lb, ub = PARAMS[topo]["lb"], PARAMS[topo]["ub"]
    if n <= 0:
        return np.empty((0, len(lb)), dtype=float)
    if RANDOM_SAMPLER in {"lhs", "latin", "latin_hypercube"} and QMC_AVAILABLE:
        sampler = qmc.LatinHypercube(d=len(lb), seed=int(rng.integers(0, 2**31 - 1)))
        u = sampler.random(n)
    else:
        u = rng.random((n, len(lb)))
    return lb + (ub - lb) * u

def normalized_param_distance(topo: str, x: np.ndarray, y: np.ndarray) -> float:
    lb, ub = PARAMS[topo]["lb"], PARAMS[topo]["ub"]
    span = np.maximum(ub - lb, EPS)
    z = (np.asarray(x, dtype=float) - np.asarray(y, dtype=float)) / span
    return float(np.linalg.norm(z) / math.sqrt(len(span)))

def select_diverse_candidates(topo: str, pool: List[Dict[str, Any]], n: int) -> List[Dict[str, Any]]:
    if not DIVERSE_REFINE or n <= 0:
        return pool[:n]
    selected: List[Dict[str, Any]] = []
    for ev in pool:
        if all(normalized_param_distance(topo, ev["x"], s["x"]) >= DIVERSE_MIN_DIST for s in selected):
            selected.append(ev)
            if len(selected) >= n:
                return selected
    for ev in pool:
        if not any(np.allclose(ev["x"], s["x"], rtol=0.0, atol=1e-9) for s in selected):
            selected.append(ev)
            if len(selected) >= n:
                break
    return selected

def elite_neighbor_samples(topo: str, pool: List[Dict[str, Any]], n_elites: int, n_per_elite: int) -> np.ndarray:
    """Create local multistart samples around good but unrefined candidates.

    This is a basin-access aid, not a statistical analysis.  It improves the
    chance that one deterministic run reaches the same low-RMS valleys that are
    otherwise found only by changing the outer random seed.
    """
    if n_elites <= 0 or n_per_elite <= 0 or not pool:
        return np.empty((0, len(PARAMS[topo]["lb"])), dtype=float)
    lb, ub = PARAMS[topo]["lb"], PARAMS[topo]["ub"]
    span = np.maximum(ub - lb, EPS)
    scales = ELITE_NEIGHBOR_SCALES or [0.08]
    elites = select_diverse_candidates(topo, pool, min(n_elites, len(pool)))
    samples: List[np.ndarray] = []
    for ev in elites:
        x0 = np.asarray(ev["x"], dtype=float)
        for j in range(n_per_elite):
            scale = scales[j % len(scales)]
            z = x0 + rng.normal(0.0, scale, size=len(x0)) * span
            samples.append(np.minimum(np.maximum(z, lb), ub))
    if not samples:
        return np.empty((0, len(lb)), dtype=float)
    return dedupe_candidate_vectors(np.asarray(samples, dtype=float))

def shape_seeded_samples(topo: str, target: np.ndarray, n: int) -> np.ndarray:
    """Deterministic, curve-scaled starting points added to the random pool.

    These are not statistics; they are a single-run basin-access aid.  The
    templates use only the target width/height to place each mechanism family
    near physically plausible link scales before DE/local refinement.
    """
    if n <= 0:
        return np.empty((0, len(PARAMS[topo]["lb"])), dtype=float)
    lb, ub = PARAMS[topo]["lb"], PARAMS[topo]["ub"]
    width = max(path_width(target), EPS)
    height = max(path_height(target), 0.25 * width)
    scale = max(width, height)
    seeds: List[np.ndarray] = []

    for k in range(n):
        a = (k + 0.5) / max(n, 1)
        jitter = 0.88 + 0.24 * rng.random()
        phi = -math.pi + 2.0 * math.pi * ((3 * k + 1) % max(n, 2)) / max(n, 2)
        if topo == "fourbar":
            x = np.array([
                (0.80 + 0.65 * a) * width,
                (0.25 + 0.45 * ((k % 5) / 4.0)) * scale,
                (0.65 + 0.80 * ((k % 7) / 6.0)) * scale,
                (0.55 + 0.70 * ((k % 6) / 5.0)) * scale,
                (0.20 + 0.95 * ((k % 8) / 7.0)) * scale,
                phi,
            ], dtype=float)
            jitter_idx = [0, 1, 2, 3, 4]
        elif topo == "slidercrank":
            x = np.array([
                (0.18 + 0.50 * ((k % 6) / 5.0)) * scale,
                (0.65 + 1.35 * a) * scale,
                -0.20 + 1.45 * ((k % 7) / 6.0),
                (-0.65 + 1.30 * ((k % 8) / 7.0)) * height,
                (-0.35 + 0.70 * ((k % 5) / 4.0)) * width,
            ], dtype=float)
            jitter_idx = [0, 1, 3, 4]
        elif topo in {"wattI", "wattII"}:
            x = np.array([
                (0.85 + 0.80 * a) * width,
                (0.25 + 0.50 * ((k % 6) / 5.0)) * scale,
                (0.65 + 0.85 * ((k % 7) / 6.0)) * scale,
                (0.55 + 0.75 * ((k % 5) / 4.0)) * scale,
                (0.20 + 0.95 * ((k % 8) / 7.0)) * scale,
                phi,
                (0.25 + 1.45 * ((k % 9) / 8.0)) * width,
                (-0.70 + 1.40 * ((k % 10) / 9.0)) * height,
                (0.45 + 0.95 * ((k % 7) / 6.0)) * scale,
                (0.45 + 0.95 * ((k % 6) / 5.0)) * scale,
                (-0.30 + 1.30 * ((k % 9) / 8.0)) * scale,
                (-0.85 + 1.70 * ((k % 10) / 9.0)) * height,
            ], dtype=float)
            jitter_idx = [0, 1, 2, 3, 4, 6, 7, 8, 9, 10, 11]
        elif topo in {"stephensonI", "stephensonII"}:
            x = np.array([
                (0.85 + 0.80 * a) * width,
                (0.25 + 0.50 * ((k % 6) / 5.0)) * scale,
                (0.65 + 0.85 * ((k % 7) / 6.0)) * scale,
                (0.55 + 0.75 * ((k % 5) / 4.0)) * scale,
                (0.45 + 0.95 * ((k % 7) / 6.0)) * scale,
                (0.45 + 0.95 * ((k % 6) / 5.0)) * scale,
                (-0.30 + 1.30 * ((k % 9) / 8.0)) * scale,
                (-0.85 + 1.70 * ((k % 10) / 9.0)) * height,
            ], dtype=float)
            jitter_idx = [0, 1, 2, 3, 4, 5, 6, 7]
        elif topo == "stephensonIII":
            x = np.array([
                (0.85 + 0.80 * a) * width,
                (0.28 + 0.48 * ((k % 6) / 5.0)) * scale,
                (0.70 + 0.75 * ((k % 7) / 6.0)) * scale,
                (0.55 + 0.70 * ((k % 5) / 4.0)) * scale,
                (0.25 + 0.85 * ((k % 8) / 7.0)) * scale,
                phi,
                (0.40 + 0.95 * ((k % 7) / 6.0)) * scale,
                (0.45 + 0.90 * ((k % 6) / 5.0)) * scale,
                (-0.25 + 1.30 * ((k % 9) / 8.0)) * scale,
                (-0.85 + 1.70 * ((k % 10) / 9.0)) * height,
            ], dtype=float)
            jitter_idx = [0, 1, 2, 3, 4, 6, 7, 8, 9]
        elif topo == "gearedFiveBar":
            ratio_choices = [-2.0, -1.5, -1.0, 0.5, 1.0, 1.5, 2.0]
            x = np.array([
                (0.75 + 0.80 * a) * width,
                (0.20 + 0.55 * ((k % 6) / 5.0)) * scale,
                (0.20 + 0.55 * ((k % 5) / 4.0)) * scale,
                (0.45 + 1.10 * ((k % 7) / 6.0)) * scale,
                (0.45 + 1.10 * ((k % 8) / 7.0)) * scale,
                ratio_choices[k % len(ratio_choices)],
                phi,
                (0.05 + 0.80 * ((k % 8) / 7.0)) * scale,
                -phi,
            ], dtype=float)
            jitter_idx = [0, 1, 2, 3, 4, 7]
        else:
            x = random_sample_topology(topo, 1)[0]
            jitter_idx = list(range(len(x)))
        x[jitter_idx] *= jitter
        seeds.append(np.minimum(np.maximum(x, lb), ub))
    return np.asarray(seeds, dtype=float)

def dedupe_candidate_vectors(X: np.ndarray, decimals: int = 5) -> np.ndarray:
    if len(X) == 0:
        return X
    _, keep = np.unique(np.round(X, decimals=decimals), axis=0, return_index=True)
    keep.sort()
    return X[keep]

def search_topology(topo: str, target: np.ndarray, theta: np.ndarray) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]]]:
    print(f"\n--- Searching {topo} ---")
    all_valid: List[Dict[str, Any]] = []
    pool: List[Dict[str, Any]] = []

    # Single-run basin access: shape-scaled templates plus random near misses.
    X_seed = shape_seeded_samples(topo, target, SHAPE_SEEDS_PER_TOPOLOGY)
    X_rand = random_sample_topology(topo, N_RANDOM_PER_TOPOLOGY)
    X = dedupe_candidate_vectors(np.vstack([X_seed, X_rand]))
    for i, x in enumerate(X):
        ev = evaluate_candidate(topo, x, target, theta)
        if ev["valid"]:
            all_valid.append(ev)
            pool.append(ev)
    pool.sort(key=lambda e: e["cost"])
    print(f"sampled={len(X)} ({len(X_seed)} shape seeds + {len(X_rand)} random) | valid={len(all_valid)} | best RMS={pool[0]['rms_mm']:.3f} mm" if pool else "no valid sampled candidates")

    # Local neighborhoods around elite unrefined samples make a single run less
    # dependent on the outer seed while keeping the search within the same
    # mechanism family.
    if pool and ELITE_NEIGHBOR_SAMPLES > 0:
        X_near = elite_neighbor_samples(topo, pool, ELITE_NEIGHBORHOODS, ELITE_NEIGHBOR_SAMPLES)
        near_valid = 0
        for x in X_near:
            ev = evaluate_candidate(topo, x, target, theta)
            if ev["valid"]:
                near_valid += 1
                all_valid.append(ev)
                pool.append(ev)
        pool.sort(key=lambda e: e["cost"])
        print(
            f"elite neighborhoods={len(X_near)} | valid={near_valid} | "
            f"best RMS={pool[0]['rms_mm']:.3f} mm"
        )

    # DE global best
    if RUN_DE_GLOBAL:
        f = objective_for_topology(topo, target, theta)
        bounds = list(zip(PARAMS[topo]["lb"], PARAMS[topo]["ub"]))
        for restart in range(max(1, DE_RESTARTS)):
            try:
                label = "" if DE_RESTARTS == 1 else f" {restart + 1}/{DE_RESTARTS}"
                print(f"running DE global{label}...")
                init = "latinhypercube" if RANDOM_SAMPLER in {"lhs", "latin", "latin_hypercube"} else "random"
                seed = int(rng.integers(0, 2**31 - 1))
                if restart > 0 and pool:
                    seed_vectors = [ev["x"] for ev in select_diverse_candidates(topo, pool, min(DE_POPSIZE, len(pool)))]
                    random_fill = random_sample_topology(topo, max(0, DE_POPSIZE * len(bounds) - len(seed_vectors)))
                    init = np.vstack([np.asarray(seed_vectors, dtype=float), random_fill])
                    if len(init) < 5:
                        init = "latinhypercube"
                de = differential_evolution(
                    f,
                    bounds=bounds,
                    maxiter=DE_MAXITER,
                    popsize=DE_POPSIZE,
                    polish=False,
                    seed=seed,
                    tol=1e-4,
                    init=init,
                    updating="immediate",
                    workers=1,
                )
                ev_de = evaluate_candidate(topo, de.x, target, theta)
                if ev_de["valid"]:
                    pool.append(ev_de)
                    all_valid.append(ev_de)
                    pool.sort(key=lambda e: e["cost"])
                    print(f"DE best RMS={ev_de['rms_mm']:.3f} MAX={ev_de['max_mm']:.3f}")
            except Exception as e:
                print(f"DE failed for {topo}: {e}")

    # Select top candidates and refine
    pool.sort(key=lambda e: e["cost"])
    selected = select_diverse_candidates(topo, pool, TOPK_RANDOM_REFINE_PER_TOPOLOGY)
    refined: List[Dict[str, Any]] = []
    for k, ev in enumerate(selected, 1):
        print(f"refining {topo} {k}/{len(selected)} | RMS={ev['rms_mm']:.3f} MAX={ev['max_mm']:.3f}")
        ref = refine_candidate(topo, ev["x"], target, theta)
        refined.append(ref)
        print(f"  -> RMS={ref['rms_mm']:.3f} MAX={ref['max_mm']:.3f} accepted={ref['accepted']}")
    return all_valid, refined

# ============================================================
# 6. Printability and STL generation
# ============================================================

def candidate_param_dict(sol: Dict[str, Any]) -> Dict[str, float]:
    topo = sol["topology"]
    return {k: float(v) for k, v in zip(PARAMS[topo]["names"], sol["x"])}

def link_lengths_for_candidate(sol: Dict[str, Any]) -> List[float]:
    d = candidate_param_dict(sol)
    topo = sol["topology"]
    if topo == "fourbar":
        return [d["L1"], d["L2"], d["L3"], d["L4"]]
    if topo == "slidercrank":
        return [d["r"], d["l"]]
    if topo in {"wattI", "wattII", "stephensonI", "stephensonII", "stephensonIII"}:
        return [d["L0"], d["L2"], d["L3"], d["L4"], d["L5"], d["L6"]]
    if topo == "gearedFiveBar":
        return [d["L0"], d["L2"], d["L4"], d["L3"], d["L5"]]
    return []

def printability_score(sol: Dict[str, Any]) -> Dict[str, Any]:
    topo = sol["topology"]
    lengths = link_lengths_for_candidate(sol)
    n_links = {
        "fourbar": 3,
        "slidercrank": 2,
        "wattI": 5,
        "wattII": 5,
        "stephensonI": 5,
        "stephensonII": 5,
        "stephensonIII": 5,
        "gearedFiveBar": 5,
    }.get(topo, 99)
    has_slider = topo == "slidercrank"
    complexity = max(0, n_links - 3) * 0.35
    slider_penalty = 0.35 if has_slider else 0.0
    min_len = min(lengths) if lengths else 0.0
    max_len = max(lengths) if lengths else 999.0
    length_penalty = 0.0
    if min_len < 14.0:
        length_penalty += 0.6
    if max_len > MAX_PRINTABLE_LINK_MM:
        length_penalty += 0.6
    not_accepted_penalty = 2.0 if not sol.get("accepted", False) else 0.0
    total = sol["rms_mm"] + 0.15 * sol["max_mm"] + complexity + slider_penalty + length_penalty + not_accepted_penalty
    return {
        "printability_score": float(total),
        "n_links": n_links,
        "has_slider": has_slider,
        "min_link_mm": float(min_len),
        "max_link_mm": float(max_len),
        "accepted": bool(sol.get("accepted", False)),
        "printable": bool(max_len <= MAX_PRINTABLE_LINK_MM and min_len >= 14.0),
    }

def ensure_stl_deps():
    if not SHAPELY_AVAILABLE:
        raise RuntimeError("shapely is required for STL generation. Install with: python -m pip install shapely")

def write_mesh_stl(vertices: np.ndarray, faces: np.ndarray, path: Path):
    if TRIMESH_AVAILABLE:
        mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=True)
        mesh.export(str(path))
        return
    # ASCII fallback
    def normal(a, b, c):
        n = np.cross(b - a, c - a)
        nn = np.linalg.norm(n)
        return n / nn if nn > 1e-12 else np.array([0.0, 0.0, 1.0])
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"solid {path.stem}\n")
        for face in faces:
            a, b, c = vertices[face[0]], vertices[face[1]], vertices[face[2]]
            n = normal(a, b, c)
            f.write(f" facet normal {n[0]} {n[1]} {n[2]}\n  outer loop\n")
            for p in (a, b, c):
                f.write(f"   vertex {p[0]} {p[1]} {p[2]}\n")
            f.write("  endloop\n endfacet\n")
        f.write(f"endsolid {path.stem}\n")

def polygon_to_extruded_mesh(poly, height: float) -> Tuple[np.ndarray, np.ndarray]:
    # Triangulate and keep triangles whose centroid lies inside the polygon.
    tris = triangulate(poly)
    verts: List[Tuple[float, float, float]] = []
    faces: List[Tuple[int, int, int]] = []
    vmap: Dict[Tuple[float, float, float], int] = {}
    def vid(x, y, z):
        key = (round(float(x), 6), round(float(y), 6), round(float(z), 6))
        if key not in vmap:
            vmap[key] = len(verts)
            verts.append(key)
        return vmap[key]
    # Top and bottom faces
    for tri in tris:
        if not poly.contains(tri.representative_point()):
            continue
        coords = list(tri.exterior.coords)[:3]
        top = [vid(x, y, height) for x, y in coords]
        bot = [vid(x, y, 0.0) for x, y in coords]
        faces.append(tuple(top))
        faces.append((bot[2], bot[1], bot[0]))
    # Side walls for exterior and holes
    def add_side_ring(coords, reverse=False):
        ring = list(coords)
        if ring[0] == ring[-1]:
            ring = ring[:-1]
        n = len(ring)
        for i in range(n):
            x1, y1 = ring[i]
            x2, y2 = ring[(i+1) % n]
            b1 = vid(x1, y1, 0.0); b2 = vid(x2, y2, 0.0)
            t1 = vid(x1, y1, height); t2 = vid(x2, y2, height)
            if not reverse:
                faces.append((b1, b2, t2)); faces.append((b1, t2, t1))
            else:
                faces.append((b1, t2, b2)); faces.append((b1, t1, t2))
    add_side_ring(poly.exterior.coords, reverse=False)
    for hole in poly.interiors:
        add_side_ring(hole.coords, reverse=True)
    return np.asarray(verts, dtype=float), np.asarray(faces, dtype=int)

def export_polygon_stl(poly, height: float, path: Path):
    verts, faces = polygon_to_extruded_mesh(poly, height)
    write_mesh_stl(verts, faces, path)

def ensure_step_deps():
    if not SHAPELY_AVAILABLE:
        raise RuntimeError("shapely is required before STEP export because the 2D profiles are built with shapely.")
    if not CADQUERY_AVAILABLE:
        raise RuntimeError(
            "CadQuery is required for editable STEP export. Install with: python -m pip install cadquery"
        )

def _polygon_geometries(poly):
    """Return a list of Polygon objects from Polygon or MultiPolygon-like geometry."""
    if poly.is_empty:
        return []
    if not getattr(poly, "is_valid", True):
        poly = poly.buffer(0)
    if poly.geom_type == "Polygon":
        return [poly]
    if poly.geom_type == "MultiPolygon":
        return list(poly.geoms)
    # Try to recover polygons from collections, if any.
    if hasattr(poly, "geoms"):
        return [g for g in poly.geoms if getattr(g, "geom_type", "") == "Polygon" and not g.is_empty]
    return []

def _clean_ring_coords(coords):
    pts = [(float(x), float(y)) for x, y in coords]
    tol = 1e-7
    if len(pts) >= 2 and math.hypot(pts[0][0] - pts[-1][0], pts[0][1] - pts[-1][1]) < tol:
        pts = pts[:-1]
    # CadQuery needs at least 3 distinct, non-degenerate points.
    cleaned = []
    for pt in pts:
        if not cleaned or math.hypot(pt[0] - cleaned[-1][0], pt[1] - cleaned[-1][1]) >= tol:
            cleaned.append(pt)
    if len(cleaned) < 3:
        return cleaned

    # Remove nearly collinear vertices; buffered Shapely profiles can contain
    # tiny segments that CadQuery/OpenCascade rejects as null geometry.
    compact = []
    for pt in cleaned:
        compact.append(pt)
        while len(compact) >= 3:
            a, b, c = compact[-3], compact[-2], compact[-1]
            ab = (b[0] - a[0], b[1] - a[1])
            bc = (c[0] - b[0], c[1] - b[1])
            cross = abs(ab[0] * bc[1] - ab[1] * bc[0])
            scale = max(math.hypot(*ab) * math.hypot(*bc), tol)
            if cross / scale < 1e-6:
                compact.pop(-2)
            else:
                break
    return compact

def _extrude_clean_polygon_points(points: List[Tuple[float, float]], height: float):
    if len(points) < 3:
        raise RuntimeError("Not enough points to extrude polygon.")
    return cq.Workplane("XY").polyline(points).close().extrude(float(height))

def polygon_to_cadquery_solid(poly, height: float):
    """Convert a shapely Polygon/MultiPolygon profile into a CadQuery extruded solid.

    This creates real B-rep solids for STEP export, which are much easier to
    edit in SolidWorks than triangulated STL meshes. Holes are cut as profiles.
    """
    ensure_step_deps()
    solids = []
    for geom in _polygon_geometries(poly):
        if geom.area <= 1e-8:
            continue
        ext = _clean_ring_coords(geom.exterior.coords)
        if len(ext) < 3:
            continue
        try:
            body = _extrude_clean_polygon_points(ext, height)
            for hole in geom.interiors:
                hpts = _clean_ring_coords(hole.coords)
                if len(hpts) >= 3:
                    cutter = _extrude_clean_polygon_points(hpts, float(height) + 2.0).translate((0, 0, -1.0))
                    body = body.cut(cutter)
            solids.append(body)
        except Exception:
            # Fallback for profiles with tiny arcs/holes that OpenCascade
            # refuses as one wire: triangulate the valid Shapely area and
            # union triangular B-rep prisms.
            for tri in triangulate(geom):
                clipped = tri.intersection(geom)
                for part in _polygon_geometries(clipped):
                    if part.area <= 1e-8:
                        continue
                    pts = _clean_ring_coords(part.exterior.coords)
                    if len(pts) >= 3:
                        solids.append(_extrude_clean_polygon_points(pts, height))
    if not solids:
        raise RuntimeError("No polygonal geometry could be converted to a STEP solid.")
    result = solids[0]
    for solid in solids[1:]:
        result = result.union(solid)
    return result

def export_polygon_step(poly, height: float, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    solid = polygon_to_cadquery_solid(poly, height)
    cq_exporters.export(solid, str(path))

def export_polygon_svg(poly, path: Path, stroke_width: float = 0.25):
    """Simple 2D SVG profile export for quick visual inspection or DXF conversion workflows."""
    ensure_stl_deps()
    geoms = _polygon_geometries(poly)
    if not geoms:
        return
    minx, miny, maxx, maxy = poly.bounds
    w, h = maxx - minx, maxy - miny
    def path_from_ring(coords):
        pts = _clean_ring_coords(coords)
        if not pts:
            return ""
        d = [f"M {pts[0][0]-minx:.4f},{maxy-pts[0][1]:.4f}"]
        for x, y in pts[1:]:
            d.append(f"L {x-minx:.4f},{maxy-y:.4f}")
        d.append("Z")
        return " ".join(d)
    paths = []
    for geom in geoms:
        paths.append(path_from_ring(geom.exterior.coords))
        for hole in geom.interiors:
            paths.append(path_from_ring(hole.coords))
    svg = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w:.4f}mm" height="{h:.4f}mm" viewBox="0 0 {w:.4f} {h:.4f}">',
        f'<g fill="none" stroke="black" stroke-width="{stroke_width}">',
    ]
    for d in paths:
        if d:
            svg.append(f'<path d="{d}"/>')
    svg.extend(['</g>', '</svg>'])
    path.write_text("\n".join(svg), encoding="utf-8")

def link_two_hole_poly(length: float, width: float, hole_diam: float):
    ensure_stl_deps()
    outer = LineString([(0, 0), (length, 0)]).buffer(width / 2.0, resolution=32, cap_style=1)
    holes = [Point(0, 0).buffer(hole_diam / 2.0, resolution=32), Point(length, 0).buffer(hole_diam / 2.0, resolution=32)]
    return outer.difference(unary_union(holes))

def ternary_plate_poly(points: List[Tuple[float, float]], width: float, hole_diam: float):
    ensure_stl_deps()
    pts = [Point(x, y).buffer(width / 2.0, resolution=24) for x, y in points]
    lines = []
    for i in range(len(points)):
        for j in range(i + 1, len(points)):
            lines.append(LineString([points[i], points[j]]).buffer(width / 2.0, resolution=16, cap_style=1))
    outer = unary_union(pts + lines)
    holes = unary_union([Point(x, y).buffer(hole_diam / 2.0, resolution=24) for x, y in points])
    return outer.difference(holes)

def base_plate_poly(points: List[Tuple[float, float]], hole_diam: float, margin: float):
    ensure_stl_deps()
    xs = [p[0] for p in points]
    ys = [p[1] for p in points]
    outer = box(min(xs)-margin, min(ys)-margin, max(xs)+margin, max(ys)+margin)
    holes = unary_union([Point(x, y).buffer(hole_diam / 2.0, resolution=24) for x, y in points])
    return outer.difference(holes)

def slider_base_poly(total_length: float, rail_width: float, hole_diam: float):
    ensure_stl_deps()
    outer = box(-BASE_MARGIN_MM, -rail_width/2, total_length + BASE_MARGIN_MM, rail_width/2)
    pivot_hole = Point(0, 0).buffer(hole_diam/2, resolution=24)
    return outer.difference(pivot_hole)

def export_best_candidate_stls(sol: Dict[str, Any], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    topo = sol["topology"]
    d = candidate_param_dict(sol)
    notes = []
    notes.append(f"Selected topology: {topo}")
    notes.append(f"RMS error: {sol['rms_mm']:.4f} mm")
    notes.append(f"MAX error: {sol['max_mm']:.4f} mm")
    notes.append(f"Parameters: {json.dumps(d, indent=2)}")
    notes.append("Recommended hardware: M3 screws or 3 mm dowel pins. Ream/adjust holes after printing.")

    if topo == "fourbar":
        L1, L2, L3, L4, rP, phiP = d["L1"], d["L2"], d["L3"], d["L4"], d["rP"], d["phiP"]
        export_polygon_stl(base_plate_poly([(0,0), (L1,0)], PIN_HOLE_DIAM_MM, BASE_MARGIN_MM), BASE_THICKNESS_MM, outdir / "base_plate_fourbar.stl")
        export_polygon_stl(link_two_hole_poly(L2, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "crank_L2.stl")
        export_polygon_stl(link_two_hole_poly(L4, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "rocker_L4.stl")
        P = (rP * math.cos(phiP), rP * math.sin(phiP))
        export_polygon_stl(ternary_plate_poly([(0,0), (L3,0), P], LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "coupler_plate_ABP.stl")
        notes.append("Assembly: ground pivots O2=(0,0), O4=(L1,0). Coupler plate holes A, B, and marker P.")

    elif topo == "slidercrank":
        r, l, q, h, x0 = d["r"], d["l"], d["q"], d["h"], d["x0"]
        total = x0 + r + l + BASE_MARGIN_MM
        export_polygon_stl(slider_base_poly(total, 22.0, PIN_HOLE_DIAM_MM), BASE_THICKNESS_MM, outdir / "base_with_slider_rail.stl")
        export_polygon_stl(link_two_hole_poly(r, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "crank_r.stl")
        P = (q * l, h)
        export_polygon_stl(ternary_plate_poly([(0,0), (l,0), P], LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "connecting_rod_with_marker.stl")
        export_polygon_stl(box(-9, -9, 9, 9).difference(Point(0,0).buffer(PIN_HOLE_DIAM_MM/2, resolution=24)), PLATE_THICKNESS_MM + 2, outdir / "slider_block.stl")
        notes.append("Assembly: use a smooth guide rail for slider. Printed rail is a concept prototype; low-friction tape or rods are recommended.")

    elif topo == "stephensonIII":
        L0, L2, L3, L4, rc, phic, L5, L6, pu, pv = [d[k] for k in PARAMS[topo]["names"]]
        export_polygon_stl(base_plate_poly([(0,0), (L0,0)], PIN_HOLE_DIAM_MM, BASE_MARGIN_MM), BASE_THICKNESS_MM, outdir / "base_plate_stephensonIII.stl")
        export_polygon_stl(link_two_hole_poly(L2, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "crank_L2.stl")
        export_polygon_stl(link_two_hole_poly(L4, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "rocker_L4.stl")
        C = (rc * math.cos(phic), rc * math.sin(phic))
        export_polygon_stl(ternary_plate_poly([(0,0), (L3,0), C], LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "ternary_coupler_ABC.stl")
        P = (pu, pv)
        export_polygon_stl(ternary_plate_poly([(0,0), (L5,0), P], LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "output_link_CDP.stl")
        export_polygon_stl(link_two_hole_poly(L6, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "connector_D_O4_L6.stl")
        notes.append("Assembly: O2 and O4 are fixed on base. Coupler plate A-B-C; output link plate C-D-P; connector D-O4.")
    else:
        notes.append(f"Detailed STL export is not implemented yet for topology {topo}.")
        notes.append("The synthesis metrics and parameter set are valid; CAD plates must be generated in a later topology-specific exporter pass.")
        (outdir / "CAD_EXPORT_NOT_IMPLEMENTED.txt").write_text("\n".join(notes), encoding="utf-8")
        return

    # Common pin boss STL
    if SHAPELY_AVAILABLE:
        boss = Point(0,0).buffer(PIN_BOSS_DIAM_MM/2, resolution=32).difference(Point(0,0).buffer(PIN_HOLE_DIAM_MM/2, resolution=24))
        export_polygon_stl(boss, PLATE_THICKNESS_MM, outdir / "optional_spacer_washer.stl")

    (outdir / "BOM_and_assembly_notes.txt").write_text("\n".join(notes), encoding="utf-8")


def export_best_candidate_steps(sol: Dict[str, Any], outdir: Path) -> None:
    """Export editable STEP CAD files for the selected printable candidate.

    STEP is the preferred format for SolidWorks because it imports as B-rep
    solids instead of triangulated meshes. This function mirrors the STL
    exporter but writes .step files through CadQuery/OpenCascade.
    """
    ensure_step_deps()
    outdir.mkdir(parents=True, exist_ok=True)
    topo = sol["topology"]
    d = candidate_param_dict(sol)
    notes = []
    notes.append(f"Selected topology: {topo}")
    notes.append(f"RMS error: {sol['rms_mm']:.4f} mm")
    notes.append(f"MAX error: {sol['max_mm']:.4f} mm")
    notes.append(f"Parameters: {json.dumps(d, indent=2)}")
    notes.append("CAD format: STEP generated with CadQuery/OpenCascade for SolidWorks editing.")
    notes.append("Recommended hardware: M3 screws or 3 mm dowel pins. Ream/adjust holes after printing.")

    if topo == "fourbar":
        L1, L2, L3, L4, rP, phiP = d["L1"], d["L2"], d["L3"], d["L4"], d["rP"], d["phiP"]
        export_polygon_step(base_plate_poly([(0,0), (L1,0)], PIN_HOLE_DIAM_MM, BASE_MARGIN_MM), BASE_THICKNESS_MM, outdir / "base_plate_fourbar.step")
        export_polygon_step(link_two_hole_poly(L2, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "crank_L2.step")
        export_polygon_step(link_two_hole_poly(L4, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "rocker_L4.step")
        P = (rP * math.cos(phiP), rP * math.sin(phiP))
        export_polygon_step(ternary_plate_poly([(0,0), (L3,0), P], LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "coupler_plate_ABP.step")
        notes.append("Assembly: ground pivots O2=(0,0), O4=(L1,0). Coupler plate holes A, B, and marker P.")

    elif topo == "slidercrank":
        r, l, q, h, x0 = d["r"], d["l"], d["q"], d["h"], d["x0"]
        total = x0 + r + l + BASE_MARGIN_MM
        export_polygon_step(slider_base_poly(total, 22.0, PIN_HOLE_DIAM_MM), BASE_THICKNESS_MM, outdir / "base_with_slider_rail.step")
        export_polygon_step(link_two_hole_poly(r, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "crank_r.step")
        P = (q * l, h)
        export_polygon_step(ternary_plate_poly([(0,0), (l,0), P], LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "connecting_rod_with_marker.step")
        export_polygon_step(box(-9, -9, 9, 9).difference(Point(0,0).buffer(PIN_HOLE_DIAM_MM/2, resolution=24)), PLATE_THICKNESS_MM + 2, outdir / "slider_block.step")
        notes.append("Assembly: use a smooth guide rail for slider. STEP base/slider are editable in SolidWorks.")

    elif topo == "stephensonIII":
        L0, L2, L3, L4, rc, phic, L5, L6, pu, pv = [d[k] for k in PARAMS[topo]["names"]]
        export_polygon_step(base_plate_poly([(0,0), (L0,0)], PIN_HOLE_DIAM_MM, BASE_MARGIN_MM), BASE_THICKNESS_MM, outdir / "base_plate_stephensonIII.step")
        export_polygon_step(link_two_hole_poly(L2, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "crank_L2.step")
        export_polygon_step(link_two_hole_poly(L4, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "rocker_L4.step")
        C = (rc * math.cos(phic), rc * math.sin(phic))
        export_polygon_step(ternary_plate_poly([(0,0), (L3,0), C], LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "ternary_coupler_ABC.step")
        P = (pu, pv)
        export_polygon_step(ternary_plate_poly([(0,0), (L5,0), P], LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "output_link_CDP.step")
        export_polygon_step(link_two_hole_poly(L6, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), PLATE_THICKNESS_MM, outdir / "connector_D_O4_L6.step")
        notes.append("Assembly: O2 and O4 are fixed on base. Coupler plate A-B-C; output link plate C-D-P; connector D-O4.")
    else:
        notes.append(f"Detailed STEP export is not implemented yet for topology {topo}.")
        notes.append("The synthesis metrics and parameter set are valid; CAD plates must be generated in a later topology-specific exporter pass.")
        (outdir / "STEP_EXPORT_NOT_IMPLEMENTED.txt").write_text("\n".join(notes), encoding="utf-8")
        return

    boss = Point(0,0).buffer(PIN_BOSS_DIAM_MM/2, resolution=32).difference(Point(0,0).buffer(PIN_HOLE_DIAM_MM/2, resolution=24))
    export_polygon_step(boss, PLATE_THICKNESS_MM, outdir / "optional_spacer_washer.step")
    (outdir / "STEP_BOM_and_assembly_notes.txt").write_text("\n".join(notes), encoding="utf-8")

def export_best_candidate_svg_profiles(sol: Dict[str, Any], outdir: Path) -> None:
    """Export 2D SVG outlines of the same printable parts.

    These are useful for quick inspection or for workflows where SolidWorks/DraftSight
    can import 2D profiles before extrusion.
    """
    ensure_stl_deps()
    outdir.mkdir(parents=True, exist_ok=True)
    topo = sol["topology"]
    d = candidate_param_dict(sol)

    if topo == "fourbar":
        L1, L2, L3, L4, rP, phiP = d["L1"], d["L2"], d["L3"], d["L4"], d["rP"], d["phiP"]
        export_polygon_svg(base_plate_poly([(0,0), (L1,0)], PIN_HOLE_DIAM_MM, BASE_MARGIN_MM), outdir / "base_plate_fourbar.svg")
        export_polygon_svg(link_two_hole_poly(L2, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), outdir / "crank_L2.svg")
        export_polygon_svg(link_two_hole_poly(L4, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), outdir / "rocker_L4.svg")
        P = (rP * math.cos(phiP), rP * math.sin(phiP))
        export_polygon_svg(ternary_plate_poly([(0,0), (L3,0), P], LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), outdir / "coupler_plate_ABP.svg")
    elif topo == "slidercrank":
        r, l, q, h, x0 = d["r"], d["l"], d["q"], d["h"], d["x0"]
        total = x0 + r + l + BASE_MARGIN_MM
        export_polygon_svg(slider_base_poly(total, 22.0, PIN_HOLE_DIAM_MM), outdir / "base_with_slider_rail.svg")
        export_polygon_svg(link_two_hole_poly(r, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), outdir / "crank_r.svg")
        P = (q * l, h)
        export_polygon_svg(ternary_plate_poly([(0,0), (l,0), P], LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), outdir / "connecting_rod_with_marker.svg")
        export_polygon_svg(box(-9, -9, 9, 9).difference(Point(0,0).buffer(PIN_HOLE_DIAM_MM/2, resolution=24)), outdir / "slider_block.svg")
    elif topo == "stephensonIII":
        L0, L2, L3, L4, rc, phic, L5, L6, pu, pv = [d[k] for k in PARAMS[topo]["names"]]
        export_polygon_svg(base_plate_poly([(0,0), (L0,0)], PIN_HOLE_DIAM_MM, BASE_MARGIN_MM), outdir / "base_plate_stephensonIII.svg")
        export_polygon_svg(link_two_hole_poly(L2, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), outdir / "crank_L2.svg")
        export_polygon_svg(link_two_hole_poly(L4, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), outdir / "rocker_L4.svg")
        C = (rc * math.cos(phic), rc * math.sin(phic))
        export_polygon_svg(ternary_plate_poly([(0,0), (L3,0), C], LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), outdir / "ternary_coupler_ABC.svg")
        P = (pu, pv)
        export_polygon_svg(ternary_plate_poly([(0,0), (L5,0), P], LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), outdir / "output_link_CDP.svg")
        export_polygon_svg(link_two_hole_poly(L6, LINK_WIDTH_MM, PIN_HOLE_DIAM_MM), outdir / "connector_D_O4_L6.svg")
    else:
        (outdir / "SVG_EXPORT_NOT_IMPLEMENTED.txt").write_text(
            f"Detailed SVG export is not implemented yet for topology {topo}.\n",
            encoding="utf-8",
        )


# ============================================================
# 7. Output and main
# ============================================================

def candidate_record(sol: Dict[str, Any]) -> Dict[str, Any]:
    rec = {
        "topology": sol["topology"],
        "valid": bool(sol.get("valid", False)),
        "accepted": bool(sol.get("accepted", False)),
        "relaxed": bool(sol.get("relaxed", False)),
        "cost": sol.get("cost", np.nan),
        "rms_mm": sol.get("rms_mm", np.nan),
        "max_mm": sol.get("max_mm", np.nan),
        "smoothness": sol.get("smoothness", np.nan),
        "min_delta": sol.get("min_delta", np.nan),
        "shift": sol.get("shift", np.nan),
        "before_rms_mm": sol.get("before_rms_mm", np.nan),
        "before_max_mm": sol.get("before_max_mm", np.nan),
        "refine_nfev": sol.get("refine_nfev", np.nan),
        "ls_polish_nfev": sol.get("ls_polish_nfev", np.nan),
    }
    rec.update(candidate_param_dict(sol))
    rec.update(printability_score(sol))
    return rec

def plot_overlay(target: np.ndarray, sol: Dict[str, Any], path: Path):
    plt.figure(figsize=(7,5))
    plt.plot(target[:,0], target[:,1], "k-", lw=2.5, label="Target")
    if sol.get("P_aligned") is not None:
        plt.plot(sol["P_aligned"][:,0], sol["P_aligned"][:,1], "--", lw=2, label=f"{sol['topology']} fit")
    plt.axis("equal")
    plt.grid(True, alpha=0.25)
    plt.xlabel("x [mm]")
    plt.ylabel("y [mm]")
    plt.legend()
    plt.title("Target vs selected mechanism path")
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()

def plot_archive(df: pd.DataFrame, path: Path):
    if df.empty:
        return
    plt.figure(figsize=(7,5))
    for topo, g in df.groupby("topology"):
        plt.scatter(g["rms_mm"], g["max_mm"], s=28, alpha=0.75, label=topo)
    plt.axvline(ACCEPT_RMS_MM, ls="--", lw=1, label="RMS tol")
    plt.axhline(ACCEPT_MAX_MM, ls=":", lw=1, label="MAX tol")
    plt.xlabel("RMS error [mm]")
    plt.ylabel("MAX error [mm]")
    plt.grid(True, alpha=0.25)
    plt.legend()
    plt.title("Candidate archive")
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()

def make_output_dir() -> Path:
    safe_tag = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in OUT_TAG)
    tag = f"_{safe_tag}" if safe_tag else ""
    stamp = f"{time.strftime('%Y%m%d_%H%M%S')}_{os.getpid()}_{time.time_ns() % 1_000_000_000:09d}"
    out = OUTPUT_ROOT / f"{OUT_PREFIX}{tag}_{stamp}"
    out.mkdir(parents=True, exist_ok=True)
    if SAVE_FIGURES:
        (out / "figures").mkdir(exist_ok=True)
    return out

def save_zip(out: Path) -> Path:
    zip_path = out.with_suffix(".zip")
    if zip_path.exists():
        zip_path.unlink()
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as z:
        for f in out.rglob("*"):
            z.write(f, arcname=f.relative_to(out.parent))
    return zip_path

def main():
    t0 = time.time()
    out = make_output_dir()
    theta = theta_samples(N_POINTS)
    target, target_meta = load_or_generate_target(theta)
    pd.DataFrame(target, columns=["x", "y"]).to_csv(out / "target_path.csv", index=False)
    (out / "target_metadata.json").write_text(json.dumps(target_meta, indent=2), encoding="utf-8")

    print("="*78)
    print("FEF multi-topology synthesis for additive-manufactured path-generator prototype")
    print("="*78)
    print(f"Application: {APPLICATION_NAME}")
    print(f"Topologies: {TOPOLOGIES}")
    print(f"Geometry profile: {GEOMETRY_PROFILE} | max printable link={MAX_PRINTABLE_LINK_MM} mm")
    print(f"Tolerance: RMS <= {ACCEPT_RMS_MM} mm | MAX <= {ACCEPT_MAX_MM} mm")
    print(f"Target width: {TARGET_WIDTH_MM} mm | source: {target_meta['target_source']}")
    print(
        f"Single-run aids: sampler={RANDOM_SAMPLER} | shape_seeds/topology={SHAPE_SEEDS_PER_TOPOLOGY} | "
        f"elite_neighbors={ELITE_NEIGHBORHOODS}x{ELITE_NEIGHBOR_SAMPLES} | "
        f"DE restarts={max(1, DE_RESTARTS)} | LS polish={RUN_LS_POLISH}"
    )
    print(f"Error metric: generated path arclength-resampled={RESAMPLE_GENERATED_FOR_ERROR}")
    print(f"C++ geometry core available: {CPP_CORE_AVAILABLE}")
    if not CPP_CORE_AVAILABLE and CPP_CORE_IMPORT_ERROR:
        print(f"C++ geometry core fallback reason: {CPP_CORE_IMPORT_ERROR}")
    print(f"Numba available: {NUMBA_AVAILABLE}")
    print(f"Shapely available: {SHAPELY_AVAILABLE} | Trimesh available: {TRIMESH_AVAILABLE} | CadQuery available: {CADQUERY_AVAILABLE}")

    all_valid: List[Dict[str, Any]] = []
    refined: List[Dict[str, Any]] = []
    for topo in TOPOLOGIES:
        valid, ref = search_topology(topo, target, theta)
        all_valid.extend(valid)
        refined.extend(ref)

    refined.sort(key=lambda e: e["rms_mm"])
    df_ref = pd.DataFrame([candidate_record(e) for e in refined])
    df_ref.to_csv(out / "candidate_archive.csv", index=False)
    df_valid = pd.DataFrame([candidate_record(e) for e in all_valid if e.get("valid", False)])
    df_valid.to_csv(out / "all_valid_evaluations.csv", index=False)

    if not refined:
        print("No refined candidates found.")
        if EXPORT_ZIP:
            zip_path = save_zip(out)
            print("ZIP:", zip_path)
        return

    best_accuracy = min(refined, key=lambda e: e["rms_mm"])
    best_printable = min(refined, key=lambda e: printability_score(e)["printability_score"])

    best_accuracy_json = {"metrics": candidate_record(best_accuracy), "params": candidate_param_dict(best_accuracy), "target_meta": target_meta}
    best_printable_json = {"metrics": candidate_record(best_printable), "params": candidate_param_dict(best_printable), "target_meta": target_meta}
    (out / "best_accuracy.json").write_text(json.dumps(best_accuracy_json, indent=2, default=float), encoding="utf-8")
    (out / "best_printable.json").write_text(json.dumps(best_printable_json, indent=2, default=float), encoding="utf-8")

    if SAVE_FIGURES:
        plot_overlay(target, best_accuracy, out / "figures" / "best_accuracy_overlay.png")
        plot_overlay(target, best_printable, out / "figures" / "best_printable_overlay.png")
        plot_archive(df_ref, out / "figures" / "archive_rms_max.png")

    # Export STLs for slicing and STEP files for SolidWorks editing.
    if EXPORT_STL:
        stl_dir = out / "stl_best_printable"
        try:
            export_best_candidate_stls(best_printable, stl_dir)
            stl_status = "STL export completed"
        except Exception as e:
            stl_status = f"STL export failed: {e}"
            (out / "STL_EXPORT_ERROR.txt").write_text(stl_status, encoding="utf-8")
    else:
        stl_status = "STL export skipped"

    if EXPORT_STEP:
        step_dir = out / "step_best_printable"
        try:
            export_best_candidate_steps(best_printable, step_dir)
            step_status = "STEP export completed"
        except Exception as e:
            step_status = f"STEP export failed: {e}"
            (out / "STEP_EXPORT_ERROR.txt").write_text(step_status, encoding="utf-8")
    else:
        step_status = "STEP export skipped"

    if EXPORT_SVG:
        svg_dir = out / "svg_profiles_best_printable"
        try:
            export_best_candidate_svg_profiles(best_printable, svg_dir)
            svg_status = "SVG profile export completed"
        except Exception as e:
            svg_status = f"SVG profile export failed: {e}"
            (out / "SVG_EXPORT_ERROR.txt").write_text(svg_status, encoding="utf-8")
    else:
        svg_status = "SVG profile export skipped"

    summary = {
        "seed": SEED,
        "application": APPLICATION_NAME,
        "target_meta": target_meta,
        "n_points": N_POINTS,
        "topologies": TOPOLOGIES,
        "geometry_profile": GEOMETRY_PROFILE,
        "max_printable_link_mm": MAX_PRINTABLE_LINK_MM,
        "n_random_per_topology": N_RANDOM_PER_TOPOLOGY,
        "shape_seeds_per_topology": SHAPE_SEEDS_PER_TOPOLOGY,
        "random_sampler": RANDOM_SAMPLER,
        "qmc_available": QMC_AVAILABLE,
        "cpp_core_available": CPP_CORE_AVAILABLE,
        "cpp_core_import_error": CPP_CORE_IMPORT_ERROR,
        "elite_neighborhoods": ELITE_NEIGHBORHOODS,
        "elite_neighbor_samples": ELITE_NEIGHBOR_SAMPLES,
        "elite_neighbor_scales": ELITE_NEIGHBOR_SCALES,
        "diverse_refine": DIVERSE_REFINE,
        "diverse_min_dist": DIVERSE_MIN_DIST,
        "resample_generated_for_error": RESAMPLE_GENERATED_FOR_ERROR,
        "save_figures": SAVE_FIGURES,
        "export_stl": EXPORT_STL,
        "export_step": EXPORT_STEP,
        "export_svg": EXPORT_SVG,
        "export_zip": EXPORT_ZIP,
        "run_ls_polish": RUN_LS_POLISH,
        "ls_max_nfev": LS_MAX_NFEV,
        "de_restarts": max(1, DE_RESTARTS),
        "de_maxiter": DE_MAXITER,
        "de_popsize": DE_POPSIZE,
        "accept_rms_mm": ACCEPT_RMS_MM,
        "accept_max_mm": ACCEPT_MAX_MM,
        "best_accuracy": candidate_record(best_accuracy),
        "best_printable": candidate_record(best_printable),
        "accepted_count": int(sum(1 for e in refined if e.get("accepted", False))),
        "refined_count": int(len(refined)),
        "valid_count": int(len(all_valid)),
        "stl_status": stl_status,
        "step_status": step_status,
        "svg_status": svg_status,
        "cadquery_available": CADQUERY_AVAILABLE,
        "elapsed_seconds": time.time() - t0,
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
    }
    (out / "run_summary.json").write_text(json.dumps(summary, indent=2, default=float), encoding="utf-8")

    print("\n" + "="*78)
    print("RESULTS")
    print("="*78)
    print("Best accuracy:")
    print(json.dumps(candidate_record(best_accuracy), indent=2, default=float))
    print("\nBest printable / STL candidate:")
    print(json.dumps(candidate_record(best_printable), indent=2, default=float))
    print("\n", stl_status)
    print(step_status)
    print(svg_status)
    print("Results folder:", out)
    if EXPORT_ZIP:
        zip_path = save_zip(out)
        print("ZIP:", zip_path)
    else:
        print("ZIP: skipped")

if __name__ == "__main__":
    main()
