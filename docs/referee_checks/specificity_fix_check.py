"""Does the proposal's membership statistic fix the specificity failure found by ref1_check.py?

ref1_check.py showed that a point-null cross-fold test of "tumour t expresses program m" fires
when t carries an unrelated program w that shares a few genes with m. In that setup kappa(m) > 0
is TRUE (the shared genes are co-expressed through w), so the failure is one of the estimand,
not of calibration. The proposal (02, section 1.2 definition 1 and section 2.4) answers with a
different estimand and test:

* kappa_tm is the mean off-diagonal latent correlation of m's genes after removing the
  activities of every OTHER dictionary program (NNLS with the loadings fixed);
* t expresses m if H0: kappa_tm <= kappa*_t is rejected, with
  kappa*_t = max(kappa_min, q95 of kappa over random gene sets matched on mean expression and
  detection in the same tumour), kappa_min = 0.05, using a cell-bootstrap standard error.

Data (as ref1_check.py): 1,000 cells, 300 genes, NB size 3, baseline mean 2 counts per gene.
Program w = genes 0-39; candidate m = 40 genes of which k are w genes (genes 40-k..39 and
100..139-k). An active program multiplies its genes' means by 2.5 (effects add when w and m
are both active on a shared gene), in an independent 30% of cells. k = 30 is added for
non-member tumours to show where the floor alone stops protecting.

Scenarios: "w only" (t does not express m: every rejection is a false membership),
"m only" and "m and w" (t expresses m: rejections measure power).

Arms:

* A, naive: ref1's statistic. Exact NB thinning into two folds (true size), Pearson residuals,
  one-sided cross-fold z-test of kappa(m) > 0 at 0.05.
* B2, w unknown: kappa(m) = mean off-diagonal of ``mpaudit.coexpr.latent_corr`` (CS-CORE
  moment estimator, size factor = total UMI) on all cells; matched-random null band (q95 of
  ``n_null`` sets from ``mpaudit.coexpr.matched_random_sets``) and floor 0.05; reject if
  (kappa - kappa*) / SE_boot > 1.645. The size factor includes the programs' own counts, so
  an active w raises it (by about 20% in the flat setting); this depth normalisation pulls
  shared and non-shared genes of m apart and shrinks kappa(m) below the raw-count value
  ``kappa_marg`` that the table prints. ``--size-factor depth`` uses the true depth instead
  (the same for every cell here), which removes that shrinkage.
* B1, w known: as B2, after conditioning on w. The dictionary is {w, m} with their TRUE
  loadings and the true baseline profile (generous: a perfect dictionary); each cell's
  activities are fitted by NNLS of its expression proportions on the dictionary; w's activity
  enters the CS-CORE mean model as a covariate (``latent_corr_given``), so kappa is a partial
  latent correlation given w. The null band is computed with the same conditioning.

Baseline settings: "flat" is exactly ref1 (every gene has mean 2). Its caveat: w genes are
then the only genes with a raised mean, so mean-matched random sets for m's shared genes are
drawn from the other w genes and the null band inherits the overlap; in a tumour expressing
both programs, sets matched to m's genes are drawn from w and the unconditioned null band
reaches kappa(w), which costs B2 its power. "varied" draws gene
baselines from 2 * lognormal(sd 0.5, mean 1), so matching pools mix w and non-w genes. In both
settings the gene universe is only 300 genes, so a random 40-gene set contains about
40 * 40 / 300 = 5 w genes; with 7,000 genes it would contain fewer than one, and the null band
here is wider than it would be in real data.

Caveats: one tumour, so no BH step; the dictionary is perfect and not learned from these
cells, so all 1,000 cells act as the held-out membership fold (the proposal's 40/30/30 split
would leave 300); no private programs, ambient RNA or depth variation; power is checked only at
an effect far above the floor; the statistic is the proposal's, not a final implementation.
"""
from __future__ import annotations

import argparse
import itertools
import time
from multiprocessing import Pool

import numpy as np
from scipy.stats import norm

from mpaudit.coexpr import latent_corr, matched_random_sets

N_CELLS = 1000
N_GENES = 300
NB_SIZE = 3.0
BASE_MEAN = 2.0
FOLD_CHANGE = 2.5
ACTIVE_FRAC = 0.3
SET_SIZE = 40
KAPPA_MIN = 0.05
ALPHA = 0.05
W_GENES = np.arange(SET_SIZE)
SCENARIOS = ("w only", "m only", "m and w")
SETTINGS = ("flat", "varied")
VARIED_SD = 0.5


def m_genes(k: int) -> np.ndarray:
    """Genes of candidate m: the last ``k`` w genes plus ``SET_SIZE - k`` genes outside w."""
    return np.concatenate([np.arange(SET_SIZE - k, SET_SIZE), np.arange(100, 100 + SET_SIZE - k)])


def program_corr() -> float:
    """Latent correlation of two genes of one program: Var(f) / (E[f^2] / size + Var(f))."""
    f = np.array([1.0, FOLD_CHANGE])
    p = np.array([1 - ACTIVE_FRAC, ACTIVE_FRAC])
    ef, ef2 = (p * f).sum(), (p * f**2).sum()
    return (ef2 - ef**2) / (ef2 / NB_SIZE + ef2 - ef**2)


def marginal_kappa(k: int) -> float:
    """kappa(m) in a "w only" tumour without conditioning: only pairs of shared genes correlate."""
    return k * (k - 1) / (SET_SIZE * (SET_SIZE - 1)) * program_corr()


def simulate(rng: np.random.Generator, setting: str, scenario: str, k: int) -> tuple[np.ndarray, np.ndarray]:
    """Counts (cells x genes) and the per-gene baseline means."""
    if setting == "flat":
        base = np.full(N_GENES, BASE_MEAN)
    else:
        base = BASE_MEAN * np.exp(rng.normal(-VARIED_SD**2 / 2, VARIED_SD, N_GENES))
    effect = np.zeros((N_CELLS, N_GENES))
    if scenario in ("w only", "m and w"):
        effect[np.ix_(rng.random(N_CELLS) < ACTIVE_FRAC, W_GENES)] += FOLD_CHANGE - 1
    if scenario in ("m only", "m and w"):
        effect[np.ix_(rng.random(N_CELLS) < ACTIVE_FRAC, m_genes(k))] += FOLD_CHANGE - 1
    lam = base[None, :] * (1 + effect)
    return rng.negative_binomial(NB_SIZE, NB_SIZE / (NB_SIZE + lam)), base


def naive_cross_fold_z(X: np.ndarray, genes: np.ndarray, rng: np.random.Generator) -> float:
    """ref1_check.py's statistic: exact NB thinning into halves, cross-fold z of kappa(m) > 0."""
    A = rng.binomial(X, rng.beta(0.5 * NB_SIZE, 0.5 * NB_SIZE, size=X.shape))
    C = X - A

    def resid(Y: np.ndarray) -> np.ndarray:
        mu = Y.mean(0, keepdims=True)
        return (Y - mu) / np.sqrt(mu + mu**2 / (0.5 * NB_SIZE))

    RA, RC = resid(A), resid(C)
    v = np.zeros(X.shape[1])
    v[genes] = 1.0
    v /= np.linalg.norm(v)
    h = (RA @ v) * (RC @ v) - (RA * RC) @ (v**2)
    return float(np.sqrt(len(h)) * h.mean() / h.std())


def nnls_small(Gram: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Exact NNLS for a few programs, vectorised over cells, by enumerating active sets.

    ``Gram`` is the (p x p) Gram matrix of the loadings and ``b`` the (cells x p) loadings-data
    products. Every active set's unconstrained solution is evaluated; the feasible one with the
    lowest objective is the NNLS solution (the optimum is the unconstrained fit on its support).
    """
    n, p = b.shape
    best = np.zeros((n, p))
    best_obj = np.zeros(n)
    for r in range(1, p + 1):
        for S in itertools.combinations(range(p), r):
            S = list(S)
            theta = np.zeros((n, p))
            theta[:, S] = np.linalg.solve(Gram[np.ix_(S, S)], b[:, S].T).T
            obj = 0.5 * np.einsum("ci,ij,cj->c", theta, Gram, theta) - (b * theta).sum(1)
            better = (theta >= 0).all(1) & (obj < best_obj)
            best[better], best_obj[better] = theta[better], obj[better]
    return best


def program_activities(Y: np.ndarray, s: np.ndarray, loadings: np.ndarray, baseline: np.ndarray) -> np.ndarray:
    """Per-cell NNLS activities of fixed loadings (programs x genes) above a fixed baseline.

    Expression proportions ``e_c = y_c / s_c`` are modelled as ``baseline + theta_c @ loadings``
    with ``theta_c >= 0``. Both the loadings and the baseline profile are given (a perfect
    dictionary): with an estimated baseline, baseline and activities trade off along the loadings
    and only the non-negativity constraint pins them.
    """
    e = Y / s[:, None]
    return nnls_small(loadings @ loadings.T, (e - baseline) @ loadings.T)


def latent_corr_given(Y: np.ndarray, s: np.ndarray, Z: np.ndarray | None = None, n_iter: int = 2) -> np.ndarray:
    """``mpaudit.coexpr.latent_corr`` with cell-specific means ``s_c (alpha_g + Z_c beta_g)``.

    ``alpha_g`` and ``beta_g`` are fitted gene by gene by least squares weighted by ``1 / s_c``
    (Poisson variance), which reduces to ``mu_g = sum(y) / sum(s)`` when ``Z`` is None, so the
    function then equals ``latent_corr``. With ``Z`` the result is a latent partial correlation
    given the covariates. The moment equations and IRLS step are those of ``latent_corr``.
    """
    Y = np.asarray(Y, dtype=np.float64)
    s = np.asarray(s, dtype=np.float64)
    D = np.ones((len(s), 1)) if Z is None else np.column_stack([np.ones(len(s)), Z - Z.mean(0)])
    coef = np.linalg.solve(D.T @ (s[:, None] * D), D.T @ Y)
    m = np.maximum(s[:, None] * (D @ coef), 1e-12)
    mu = Y.sum(axis=0) / max(s.sum(), 1e-12)
    R = Y - m
    s2 = s**2
    Dm = R**2 - m
    var = (Dm * s2[:, None]).sum(axis=0) / float((s2**2).sum())
    cov = (R * s2[:, None]).T @ R / float((s2**2).sum())
    for _ in range(n_iter):
        sig = np.clip(var, 1e-12 * np.maximum(mu, 1e-12) ** 2, None)
        V = np.maximum(m + s2[:, None] * sig[None, :], 1e-12)
        var = ((Dm * s2[:, None]) / V**2).sum(axis=0) / ((s2**2)[:, None] / V**2).sum(axis=0)
        Rv = R / V
        cov = (Rv * s2[:, None]).T @ Rv / (((s2**2)[:, None] / V).T @ (1.0 / V))
    np.fill_diagonal(cov, var)
    with np.errstate(invalid="ignore", divide="ignore"):
        sd = np.sqrt(np.where(var > 0, var, np.nan))
        C = cov / np.outer(sd, sd)
    return np.clip(C, -1.0, 1.0)


def set_kappa(C: np.ndarray, genes: np.ndarray) -> float:
    """Mean off-diagonal entry of the correlation submatrix of ``genes``."""
    sub = C[np.ix_(genes, genes)]
    return float(np.nanmean(sub[np.triu_indices(len(genes), k=1)]))


def membership_test(
    Y: np.ndarray,
    s: np.ndarray,
    genes: np.ndarray,
    Z: np.ndarray | None,
    null_sets: list[np.ndarray],
    rng: np.random.Generator,
    n_boot: int,
) -> dict:
    """One-sided test of H0: kappa <= kappa* = max(KAPPA_MIN, null q95), bootstrap SE.

    The latent correlation of a set equals the submatrix of the all-gene matrix (every weight in
    the estimator is gene-wise), so kappa and the whole null band come from one matrix. The
    bootstrap is skipped when kappa <= kappa*, where H0 cannot be rejected whatever the SE.
    """
    C = latent_corr(Y, s) if Z is None else latent_corr_given(Y, s, Z)
    kap = set_kappa(C, genes)
    q95 = float(np.quantile([set_kappa(C, g) for g in null_sets], 0.95))
    kstar = max(KAPPA_MIN, q95)
    reject = False
    if kap > kstar:
        boots = []
        for _ in range(n_boot):
            idx = rng.integers(0, len(s), len(s))
            Cb = latent_corr_given(Y[np.ix_(idx, genes)], s[idx], None if Z is None else Z[idx])
            boots.append(set_kappa(Cb, np.arange(len(genes))))
        reject = (kap - kstar) / np.std(boots, ddof=1) > norm.ppf(1 - ALPHA)
    return {"kappa": kap, "q95": q95, "floor_binds": q95 < KAPPA_MIN, "reject": bool(reject)}


def one_replicate(task: tuple[int, str, str, int, int, int, int, str]) -> dict:
    """All arms on one simulated tumour."""
    seed, setting, scenario, k, rep, n_null, n_boot, size_factor = task
    rng = np.random.default_rng([seed, SETTINGS.index(setting), SCENARIOS.index(scenario), k, rep])
    X, base = simulate(rng, setting, scenario, k)
    genes = m_genes(k)
    out = {"setting": setting, "scenario": scenario, "k": k,
           "A": naive_cross_fold_z(X, genes, rng) > norm.ppf(1 - ALPHA)}
    Y = X.astype(np.float64)
    s = Y.sum(1) if size_factor == "total" else np.full(N_CELLS, base.sum())
    null_sets = matched_random_sets(genes, Y.sum(0) / s.sum(), (X > 0).mean(0), n_null, rng)
    out.update({f"B2_{key}": v for key, v in membership_test(Y, s, genes, None, null_sets, rng, n_boot).items()})
    loadings = np.zeros((2, N_GENES))
    loadings[0, W_GENES] = base[W_GENES]
    loadings[1, genes] = base[genes]
    loadings /= loadings.sum(1, keepdims=True)
    a_w = program_activities(Y, s, loadings, base / base.sum())[:, :1]
    out.update({f"B1_{key}": v for key, v in membership_test(Y, s, genes, a_w, null_sets, rng, n_boot).items()})
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--k-nonmember", type=int, nargs="+", default=[0, 2, 5, 10, 20, 30],
                    help="shared genes k for the 'w only' scenario")
    ap.add_argument("--k-member", type=int, nargs="+", default=[0, 2, 5, 10, 20],
                    help="shared genes k for the 'm only' and 'm and w' scenarios")
    ap.add_argument("--settings", nargs="+", default=list(SETTINGS), choices=SETTINGS)
    ap.add_argument("--scenarios", nargs="+", default=list(SCENARIOS), choices=SCENARIOS)
    ap.add_argument("--size-factor", default="total", choices=("total", "depth"),
                    help="B arms: total UMI per cell (as mpaudit) or the true sequencing depth, "
                         "which is the same for every cell in this toy")
    ap.add_argument("--n-null", type=int, default=200)
    ap.add_argument("--n-boot", type=int, default=50)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--seed", type=int, default=1)
    args = ap.parse_args()

    t0 = time.time()
    ks = {sc: args.k_nonmember if sc == "w only" else args.k_member for sc in SCENARIOS}
    tasks = [(args.seed, st, sc, k, r, args.n_null, args.n_boot, args.size_factor)
             for st in args.settings for sc in args.scenarios for k in ks[sc] for r in range(args.reps)]
    with Pool(args.workers) as pool:
        rows = pool.map(one_replicate, tasks, chunksize=10)

    print(f"specificity_fix_check.py: {N_CELLS} cells, {N_GENES} genes, NB size {NB_SIZE:g}, "
          f"fold change {FOLD_CHANGE:g} in {ACTIVE_FRAC:.0%} of cells, |m| = |w| = {SET_SIZE}")
    print(f"reps={args.reps} per row, n_null={args.n_null}, n_boot={args.n_boot}, "
          f"kappa_min={KAPPA_MIN}, alpha={ALPHA}, seed={args.seed}, size factor={args.size_factor}")
    print(f"latent correlation of two genes of one active program (theory): {program_corr():.3f}")
    print("columns: kappa_marg = kappa(m) of raw latent expression in a 'w only' tumour (theory, no "
          "depth normalisation); A, B2, B1 = rejection rate (false membership under 'w only', power "
          "otherwise); kap_*, q95_* = mean kappa estimate and mean null q95 per arm; flr_* = share of "
          "replicates where the floor 0.05 exceeds the null q95")
    hdr = (f"{'setting':7s} {'scenario':8s} {'k':>2s} {'kappa_marg':>10s} {'A':>5s} {'B2':>5s} {'B1':>5s}"
           f" {'kap_B2':>7s} {'q95_B2':>7s} {'flr_B2':>6s} {'kap_B1':>7s} {'q95_B1':>7s} {'flr_B1':>6s}")
    print(hdr)
    for st in args.settings:
        for sc in args.scenarios:
            for k in ks[sc]:
                g = [r for r in rows if r["setting"] == st and r["scenario"] == sc and r["k"] == k]

                def mean(key: str) -> float:
                    return float(np.mean([r[key] for r in g]))

                kmarg = f"{marginal_kappa(k):10.3f}" if sc == "w only" else f"{'':10s}"
                print(f"{st:7s} {sc:8s} {k:2d} {kmarg} {mean('A'):5.2f} {mean('B2_reject'):5.2f} "
                      f"{mean('B1_reject'):5.2f} {mean('B2_kappa'):+7.3f} {mean('B2_q95'):7.3f} "
                      f"{mean('B2_floor_binds'):6.2f} {mean('B1_kappa'):+7.3f} {mean('B1_q95'):7.3f} "
                      f"{mean('B1_floor_binds'):6.2f}", flush=True)
    print(f"Monte Carlo SE of a rate p with {args.reps} reps: sqrt(p(1-p)/{args.reps}) "
          f"(0.015 at p=0.05, 0.035 at p=0.5). Runtime {time.time() - t0:.0f} s with {args.workers} workers.")


if __name__ == "__main__":
    main()
