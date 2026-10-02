"""Referee-2 toy check of InvarMP's thinned split-LR e-value for loading invariance (D2).

Procedure as specified in the concept:
  - NB data thinning X = X1 + X2 (true phi known; generous to the method)
  - fit UNRESTRICTED model on X1 (program 1 gets separate loadings in env a)
  - on X2: numerator = X1 fit (loadings + per-cell usage) evaluated on X2 env-a cells
           denominator = restricted model (shared loadings from X1), per-cell usage
                         re-optimised on X2 by convex Poisson (KL) NNLS
  - log E = loglik_X2(num) - loglik_X2(den)  [Poisson working lik, and NB true lik]
Variants:
  'spec'    : as specified (denominator re-fits per-cell usage on X2)
  'plugin'  : denominator also uses X1 usage (NOT a valid UI e-value; shows incidental-param cost)
Also: one shared kappa error (gene-wise capture misestimate in env a, common to all programs)
under H0 for program-specific delta.
"""
import numpy as np
from scipy.special import gammaln
import sys

rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 7)
G, K = 300, 3
N0, NA = 3000, 800
PHI = 3.0
EPS = 0.5


def make_W():
    W = np.full((K, G), 0.02)
    for k in range(K):
        idx = np.arange(k * 40, k * 40 + 40)
        W[k, idx] += rng.gamma(2.0, 0.5, size=40)
    # shared highly-expressed background genes
    W[:, 200:230] += 0.3
    return W / W.sum(1, keepdims=True)


def nb_sample(mu, phi):
    return rng.negative_binomial(phi, phi / (phi + mu))


def nb_thin(X, phi, eps):
    p = rng.beta(eps * phi, (1 - eps) * phi, size=X.shape)
    X1 = rng.binomial(X, p)
    return X1, X - X1


def pois_ll(X, mu):
    mu = np.clip(mu, 1e-12, None)
    return (X * np.log(mu) - mu - gammaln(X + 1)).sum()


def nb_ll(X, mu, phi):
    mu = np.clip(mu, 1e-12, None)
    return (gammaln(X + phi) - gammaln(phi) - gammaln(X + 1)
            + X * np.log(mu / (mu + phi)) + phi * np.log(phi / (mu + phi))).sum()


def kl_fit(X, Th, Wf, mask, iters=300, fixW=False):
    """KL-NMF multiplicative updates. Wf: (R x G) loadings rows; mask: (N x R) allowed rows."""
    Th = Th * mask
    for _ in range(iters):
        M = Th @ Wf + 1e-12
        R = X / M
        Th *= (R @ Wf.T) / (Wf.sum(1)[None, :] + 1e-12)
        Th *= mask
        if not fixW:
            M = Th @ Wf + 1e-12
            R = X / M
            Wf *= (Th.T @ R) / (Th.sum(0)[:, None] + 1e-12)
    return Th, Wf


def one_rep(delta_size, n_delta_genes=20, kappa_err=0.0, umi=800):
    W = make_W()
    th = rng.gamma(1.0, 1.0, size=(N0 + NA, K))
    lib = rng.lognormal(np.log(umi), 0.4, size=N0 + NA)
    th = th / th.sum(1, keepdims=True) * lib[:, None]
    env_a = np.zeros(N0 + NA, bool)
    env_a[N0:] = True
    Wa = W.copy()
    if delta_size != 0:
        g = np.arange(0, n_delta_genes)
        Wa[0, g] *= np.exp(delta_size)  # program-specific deviation in env a
    kap = np.ones(G)
    if kappa_err > 0:
        kap = np.exp(rng.normal(0, kappa_err, G))  # shared capture error, all programs
    mu = np.where(env_a[:, None], (th @ Wa) * kap[None, :], th @ W)
    X = nb_sample(mu, PHI)
    X1, X2 = nb_thin(X, PHI, EPS)
    # unrestricted fit on X1: rows [W1, W2, W3, W1a]
    Wf = np.vstack([W, W[0:1]]) * np.exp(rng.normal(0, 0.1, (K + 1, G)))
    mask = np.zeros((N0 + NA, K + 1))
    mask[~env_a, :3] = 1
    mask[env_a, 1:] = 1
    Th0 = np.where(mask > 0, EPS * lib[:, None] / K, 0.0)
    Th1, Wf1 = kl_fit(X1.astype(float), Th0.copy(), Wf.copy(), mask, iters=400)
    A = env_a
    Xa2 = X2[A].astype(float)
    sc = (1 - EPS) / EPS
    mu_num = (Th1[A] @ Wf1) * sc
    # restricted (delta=0) denominator: loadings W1,W2,W3 from X1 fit, usage refit on X2
    Wr = Wf1[:3].copy()
    th_init = np.column_stack([Th1[A][:, 3], Th1[A][:, 1], Th1[A][:, 2]]) * sc
    th_den, _ = kl_fit(Xa2, th_init.copy(), Wr, np.ones_like(th_init), iters=400, fixW=True)
    mu_den = th_den @ Wr
    mu_plug = th_init @ Wr
    out = {}
    out['spec_pois'] = pois_ll(Xa2, mu_num) - pois_ll(Xa2, mu_den)
    out['spec_nb'] = nb_ll(Xa2, mu_num, (1 - EPS) * PHI) - nb_ll(Xa2, mu_den, (1 - EPS) * PHI)
    out['plug_pois'] = pois_ll(Xa2, mu_num) - pois_ll(Xa2, mu_plug)
    out['plug_nb'] = nb_ll(Xa2, mu_num, (1 - EPS) * PHI) - nb_ll(Xa2, mu_plug, (1 - EPS) * PHI)
    return out


def summarise(name, rows):
    keys = rows[0].keys()
    for k in keys:
        v = np.array([r[k] for r in rows])
        print(f"{name:28s} {k:10s} mean logE={v.mean():+9.1f}  sd={v.std():7.1f}  "
              f"P(E>=20)={np.mean(v >= np.log(20)):.2f}  P(E>=160)={np.mean(v >= np.log(160)):.2f}",
              flush=True)


if __name__ == "__main__":
    R = 20
    for label, d, ke in [("H0 delta=0", 0.0, 0.0),
                         ("H1 delta=log1.65 x20g", 0.5, 0.0),
                         ("H1 delta=log2.72 x20g", 1.0, 0.0),
                         ("H0 delta=0, kappa sd=0.1", 0.0, 0.10),
                         ("H0 delta=0, kappa sd=0.2", 0.0, 0.20)]:
        rows = [one_rep(d, kappa_err=ke) for _ in range(R)]
        summarise(label, rows)
