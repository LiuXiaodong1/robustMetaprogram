"""Shared helpers of the evaluation harness (evaluate_c.py, evaluate_a.py).

Nothing here imports mpaudit at module load: :func:`load_mpaudit` puts the candidate package first on
``sys.path`` and checks that every imported ``mpaudit`` module comes from it.
"""
from __future__ import annotations

import os

# one BLAS/OpenMP thread per worker process; must be set before numpy is imported
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_v] = "1"

import importlib
import json
import math
import multiprocessing as mp
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np


def resolve_pkg(pkg_dir: str) -> tuple[Path, Path]:
    """``(parent, package_dir)`` for PKG_DIR given either as the ``mpaudit`` directory itself or as a
    directory that contains ``mpaudit/``."""
    p = Path(pkg_dir).resolve()
    if (p / "__init__.py").is_file() and p.name == "mpaudit":
        return p.parent, p
    if (p / "mpaudit" / "__init__.py").is_file():
        return p, p / "mpaudit"
    raise SystemExit(f"{pkg_dir}: neither an mpaudit package directory nor a directory containing mpaudit/")


def load_mpaudit(pkg_dir: str):
    """Import ``mpaudit`` from PKG_DIR (and nowhere else). Returns ``(parent, mpaudit module)``.

    The venv has an editable install of the repository's ``src``; the candidate's parent is put at
    ``sys.path[0]`` (ahead of site-packages .pth entries) and every ``mpaudit*`` module already
    imported is dropped first. Each submodule imported later is checked by :func:`check_modules`.
    """
    parent, pkg = resolve_pkg(pkg_dir)
    for name in [m for m in sys.modules if m == "mpaudit" or m.startswith("mpaudit.")]:
        del sys.modules[name]
    sys.path[:] = [x for x in sys.path if Path(x or ".").resolve() != parent]
    sys.path.insert(0, str(parent))
    importlib.invalidate_caches()
    mpa = importlib.import_module("mpaudit")
    if Path(mpa.__file__).resolve().parent != pkg:
        raise SystemExit(f"imported mpaudit from {mpa.__file__}, not from {pkg}")
    return parent, mpa


def check_modules(pkg_dir: str) -> dict:
    """Assert every imported ``mpaudit*`` module lives under PKG_DIR; return {module: file}."""
    _, pkg = resolve_pkg(pkg_dir)
    out = {}
    for name, m in list(sys.modules.items()):
        if name == "mpaudit" or name.startswith("mpaudit."):
            f = getattr(m, "__file__", None)
            if f is None:
                continue
            if pkg not in Path(f).resolve().parents:
                raise SystemExit(f"module {name} imported from {f}, outside {pkg}")
            out[name] = str(Path(f).resolve())
    return out


# --------------------------------------------------------------------------- statistics

def ratio_se(num: np.ndarray, den: np.ndarray) -> tuple[float, float, int, int]:
    """Ratio estimator ``p = sum(num) / sum(den)`` over clusters, with a cluster-robust (linearised)
    Monte Carlo standard error. Returns ``(p, se, n = sum(den), n_clusters with den > 0)``."""
    num = np.asarray(num, float)
    den = np.asarray(den, float)
    keep = den > 0
    num, den = num[keep], den[keep]
    n, k = float(den.sum()), int(keep.sum())
    if n == 0:
        return float("nan"), float("nan"), 0, 0
    p = float(num.sum() / n)
    if k < 2:
        return p, float("nan"), int(n), k
    r = num - p * den
    se = math.sqrt(k / (k - 1) * float((r**2).sum())) / n
    return p, se, int(n), k


def binom(k: int, n: int) -> dict:
    """Proportion with its binomial Monte Carlo SE (independent trials)."""
    if n == 0:
        return {"p": None, "se": None, "n": 0, "k": 0}
    p = k / n
    return {"p": round(p, 5), "se": round(math.sqrt(p * (1 - p) / n), 5), "n": int(n), "k": int(k)}


def rate(num, den) -> dict:
    p, se, n, k = ratio_se(num, den)
    return {"p": None if not np.isfinite(p) else round(p, 5), "se": None if not np.isfinite(se) else round(se, 5),
            "n": n, "n_clusters": k, "k": int(round(float(np.sum(np.asarray(num, float)[np.asarray(den, float) > 0]))))}


def p_at_least(m: int, M: int, p: float) -> float:
    """P(Binomial(M, p) >= m)."""
    from scipy import stats
    return float(stats.binom.sf(m - 1, M, p))


# --------------------------------------------------------------------------- parallel map

def _run_job(fn_name: str, module: str, pkg_dir: str, job: dict):
    t0 = time.time()
    try:
        mod = sys.modules[module] if module in sys.modules else importlib.import_module(module)
        res = getattr(mod, fn_name)(pkg_dir, job)
        return {"job": job, "ok": True, "result": res, "seconds": time.time() - t0}
    except Exception:  # keep the run going; the error is reported in the output
        return {"job": job, "ok": False, "error": traceback.format_exc(), "seconds": time.time() - t0}


def pool_map(fn_name: str, module: str, pkg_dir: str, jobs: list[dict], n_workers: int, log=print) -> list[dict]:
    """Run ``module.fn_name(pkg_dir, job)`` for every job on a fork pool, longest (``job['cost']``) first.

    The results come back in the order of ``jobs`` (not in the order the jobs finish), so every list built
    from them is the same in every run."""
    order = sorted(range(len(jobs)), key=lambda i: -jobs[i].get("cost", 1.0))  # stable: ties keep the jobs' order
    out: list[dict | None] = [None] * len(jobs)
    if n_workers <= 1:
        for i in order:
            out[i] = r = _run_job(fn_name, module, pkg_dir, jobs[i])
            log(f"  done {jobs[i].get('name')} in {r['seconds']:.0f}s{'' if r['ok'] else ' (ERROR)'}")
        return out
    ctx = mp.get_context("fork")
    with ProcessPoolExecutor(max_workers=n_workers, mp_context=ctx) as ex:
        futs = {ex.submit(_run_job, fn_name, module, pkg_dir, jobs[i]): i for i in order}
        for f in as_completed(futs):
            out[futs[f]] = r = f.result()
            log(f"  done {r['job'].get('name')} in {r['seconds']:.0f}s{'' if r['ok'] else ' (ERROR)'}", flush=True)
    return out


def dump(obj, path: str):
    def default(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return None if not np.isfinite(o) else float(o)
        if isinstance(o, (np.bool_,)):
            return bool(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        return str(o)

    def clean(o):
        if isinstance(o, float) and not math.isfinite(o):
            return None
        if isinstance(o, dict):
            return {str(k): clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        return o

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(clean(obj), indent=1, default=default))
