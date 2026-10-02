"""Referee-2 checks of CALIPER's calibration claims (toy simulations, not publishable).

Sim 1: type I of the cross-fold (CF) test with a selection-heavy engine when the null
       (no cross-gene co-expression) holds but per-gene marginals are NOT negative binomial.
Sim 2: specificity of tumour membership: a tumour expressing program P but not m,
       where m shares a few genes with P.
"""
import numpy as np
from scipy.stats import norm
import sys

rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
EPS = np.array([0.5, 0.2, 0.3])


def dirmult_thin(X, btil):
    """Split X (n x G) into 3 folds by Dirichlet-multinomial with alpha_k = eps_k * btil_j."""
    n, G = X.shape
    alpha = EPS[None, None, :] * btil[None, :, None]  # 1 x G x 3
    g = rng.gamma(np.broadcast_to(alpha, (n, G, 3)))
    p = g / g.sum(axis=2, keepdims=True)
    p = np.nan_to_num(p, nan=1 / 3)
    XA = rng.binomial(X, p[:, :, 0])
    rem = X - XA
    pb = p[:, :, 1] / np.clip(p[:, :, 1] + p[:, :, 2], 1e-12, None)
    XB = rng.binomial(rem, np.clip(pb, 0, 1))
    XC = rem - XB
    return XA, XB, XC


def resid(Xk, l, btil_k):
    n = Xk.shape[0]
    a = Xk.sum(0) / l.sum()
    nu = l[:, None] * a[None, :]
    v = nu + nu ** 2 / btil_k[None, :]
    r = (Xk - nu) / np.sqrt(np.clip(v, 1e-12, None))
    return np.clip(r, -np.sqrt(n), np.sqrt(n))


def b_mom(X, l):
    a = X.sum(0) / l.sum()
    nu = l[:, None] * a[None, :]
    den = ((X - nu) ** 2 - nu).sum(0)
    num = (nu ** 2).sum(0)
    b = np.where(den > 0, num / np.clip(den, 1e-12, None), 1e4)
    return np.clip(b, 0.05, 1e4)


def engine_pca_top(RA, k=30):
    Rc = RA - RA.mean(0)
    U, S, Vt = np.linalg.svd(Rc, full_matrices=False)
    v = Vt[0]
    if v.sum() < 0:
        v = -v
    idx = np.argsort(-v)[:k]
    m = np.zeros(RA.shape[1])
    m[idx] = np.clip(v[idx], 1e-6, None)
    return m / np.linalg.norm(m)


def engine_seed_coexp(RA, k=30):
    """More tail-driven selection: seed = gene with largest fold-A variance; take k-1 genes most
    correlated with it on fold A."""
    var = RA.var(0)
    s = np.argmax(var)
    c = (RA - RA.mean(0)).T @ (RA[:, s] - RA[:, s].mean())
    c[s] = np.inf
    idx = np.argsort(-c)[:k]
    m = np.zeros(RA.shape[1])
    m[idx] = 1.0
    return m / np.linalg.norm(m)


def cf_z(RA, RC, m):
    sA = RA @ m
    sC = RC @ m
    h = sA * sC - (RA * RC) @ (m ** 2)
    n = len(h)
    return np.sqrt(n) * h.mean() / h.std(ddof=1)


def simulate_null(scen, n=300, G=400):
    l = rng.lognormal(0, 0.4, n)
    l = l / l.mean()
    a = np.exp(rng.normal(np.log(1.0), 1.0, G))
    a = np.clip(a, 0.05, 30)
    base = l[:, None] * a[None, :]
    if scen == "nb_exact":
        b = 3.0
        mu = base
        X = rng.negative_binomial(b, b / (b + mu))
        btrue = np.full(G, b)
    elif scen == "onoff":
        # per-gene independent on/off latent state (no cross-gene structure => H0 true)
        p_on, f = 0.08, 8.0
        z = rng.random((n, G)) < p_on
        mu = base * (1 + (f - 1) * z) / (1 + (f - 1) * p_on)
        b = 10.0
        X = rng.negative_binomial(b, b / (b + mu))
        btrue = None
    elif scen == "lognormal":
        s = 1.0
        mu = base * np.exp(rng.normal(-s * s / 2, s, (n, G)))
        X = rng.poisson(mu)
        btrue = None
    return X, l, btrue


def run_sim1(scen, mode, engine, reps):
    rej = 0
    zs = []
    for _ in range(reps):
        X, l, btrue = simulate_null(scen)
        if mode == "exact":
            bt = btrue
        elif mode == "mom":
            bt = b_mom(X, l)
        elif mode == "mom_half":
            bt = 0.5 * b_mom(X, l)
        elif mode == "poisson":
            bt = np.full(X.shape[1], 1e6)
        XA, XB, XC = dirmult_thin(X, bt)
        RA = resid(XA, l, EPS[0] * bt)
        RC = resid(XC, l, EPS[2] * bt)
        m = engine(RA)
        z = cf_z(RA, RC, m)
        zs.append(z)
        rej += z > norm.ppf(0.95)
    zs = np.array(zs)
    return rej / reps, zs.mean(), zs.std()


def run_sim2(n, n_shared, reps, k=40, fc=2.0, frac=0.3, G=400):
    """Tumour expresses P (genes 0..k-1). Test fixed m = n_shared genes of P + (k-n_shared) non-P genes.
    Exact NB thinning (b known), so any rejection is 'kappa(m)>0' being true via overlap."""
    rej = 0
    rej_emp = 0
    b = 3.0
    for _ in range(reps):
        l = rng.lognormal(0, 0.4, n)
        l = l / l.mean()
        a = np.clip(np.exp(rng.normal(0, 1.0, G)), 0.05, 30)
        act = rng.random(n) < frac
        mult = np.ones((n, G))
        mult[np.ix_(act, np.arange(k))] = fc
        mu = l[:, None] * a[None, :] * mult
        X = rng.negative_binomial(b, b / (b + mu))
        bt = np.full(G, b)
        XA, XB, XC = dirmult_thin(X, bt)
        RA = resid(XA, l, EPS[0] * bt)
        RC = resid(XC, l, EPS[2] * bt)
        m = np.zeros(G)
        m[:n_shared] = 1
        m[k:k + (k - n_shared)] = 1
        m /= np.linalg.norm(m)
        z = cf_z(RA, RC, m)
        rej += z > norm.ppf(0.95)
        # empirical null: 200 pseudo-programs of k random non-P genes
        nullz = []
        nonP = np.arange(k, G)
        for _b in range(200):
            mm = np.zeros(G)
            mm[rng.choice(nonP, k, replace=False)] = 1
            mm /= np.linalg.norm(mm)
            nullz.append(cf_z(RA, RC, mm))
        nullz = np.array(nullz)
        zstar = (z - nullz.mean()) / nullz.std()
        rej_emp += zstar > norm.ppf(0.95)
    return rej / reps, rej_emp / reps


if __name__ == "__main__":
    which = sys.argv[1]
    if which == "sim1":
        reps = 300
        for scen, modes in [("nb_exact", ["exact", "mom", "poisson"]),
                            ("onoff", ["mom", "mom_half", "poisson"]),
                            ("lognormal", ["mom", "mom_half", "poisson"])]:
            for eng_name, eng in [("pca_top30", engine_pca_top), ("seed_coexp30", engine_seed_coexp)]:
                for mode in modes:
                    r = reps if mode != "poisson" else 100
                    t1, mz, sz = run_sim1(scen, mode, eng, r)
                    print(f"{scen:10s} {eng_name:13s} {mode:9s} typeI={t1:.3f} meanZ={mz:+.2f} sdZ={sz:.2f} reps={r}", flush=True)
    elif which == "sim2":
        for n in [300, 1000, 3000]:
            for ns in [0, 2, 4, 8]:
                r1, r2 = run_sim2(n, ns, reps=100 if n < 3000 else 50)
                print(f"n={n:5d} shared={ns}/40  reject(kappa>0)={r1:.2f}  reject_after_empirical_null={r2:.2f}", flush=True)
