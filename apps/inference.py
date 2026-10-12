"""Inference for a pooled combination of J estimators computed on one sample.

Input throughout: `psi` (J estimates) and `Phi` (N x J matrix of per-unit influence
contributions, scaled so that psi_j - target_j ~ sum_i Phi[i, j]). Then
    Sigma = Phi' Phi            estimated covariance matrix of the J estimators,
    v     = diag(Sigma)         their variances,
and every variance of a weighted combination w'psi is w' Sigma w. This is the
observation-level sandwich of the paper (Section 5), written without choosing a
global n: with N i.i.d. units, Sigma = Omega-hat / N.

Nothing here refits a model; it only post-processes (psi, Phi).
"""
import numpy as np

Z = {0.90: 1.6448536269514722, 0.95: 1.959963984540054}


# ---------------------------------------------------------------- dispersion
def tau2_pairwise(psi, v):
    """Pairwise-difference estimator, eq. (tau2-pairwise) of the paper."""
    i, j = np.triu_indices(len(psi), 1)
    return max(0.0, float(np.mean(((psi[i] - psi[j]) ** 2 - (v[i] + v[j])) / 2)))


def tau2_dl(psi, v):
    """DerSimonian-Laird moment estimator (reported for comparison)."""
    w = 1 / v
    pbar = np.sum(w * psi) / np.sum(w)
    q = np.sum((psi - pbar) ** 2 / v)
    c = np.sum(w) - np.sum(w ** 2) / np.sum(w)
    return max(0.0, float((q - (len(psi) - 1)) / c))


def tau2_thresholded(tau2, v, n_units, kappa=None):
    """Proposition 3(iii) on a scale-free footing: zero unless tau2 > kappa_n * mean(v).

    The paper writes the rule as n * tau2 > kappa_n, which depends on the unit of the
    outcome. Since v_j = omega_jj / n, comparing tau2 with kappa_n * mean(v) is the same
    rule with kappa_n multiplied by the average of the omega_jj: it satisfies the same
    two conditions (kappa_n -> infinity, kappa_n / n -> 0) and is invariant to rescaling.
    Default kappa_n = log n.
    """
    kappa = np.log(n_units) if kappa is None else kappa
    ratio = tau2 / float(np.mean(v))
    return (tau2 if ratio > kappa else 0.0), bool(ratio > kappa), ratio


def dispersion_z(psi, Sigma):
    """Standardised dispersion of the estimates under exact identification.

    S2 is the sample variance of the J estimates (the pairwise estimator is max(0, S2 - mean v)).
    If all targets coincide, S2 = z'Cz/(J-1) with z ~ N(0, Sigma) and C the centring matrix,
    so E(S2) = tr(C Sigma)/(J-1) and Var(S2) = 2 tr((C Sigma)^2)/(J-1)^2. With positively
    correlated estimators E(S2) is below mean v, which is why the pairwise estimator is so
    often zero. Returned: (S2 - E S2) / sd(S2), a diagnostic that uses the covariances.
    """
    J = len(psi)
    C = np.eye(J) - 1.0 / J
    M = C @ Sigma
    s2 = float(psi @ C @ psi) / (J - 1)
    return (s2 - np.trace(M) / (J - 1)) / (np.sqrt(2 * np.trace(M @ M)) / (J - 1))


# ---------------------------------------------------------------- weights
def eb_weights(v, tau2):
    w = 1 / (v + tau2)
    return w / w.sum()


def gls_weights(Sigma, rcond=1e-8):
    """Minimum-variance weights Sigma^+ 1 / 1' Sigma^+ 1 (pseudo-inverse; reported with its rank)."""
    ev = np.linalg.eigvalsh(Sigma)
    rank = int((ev > rcond * ev.max()).sum())
    g = np.linalg.pinv(Sigma, rcond=rcond, hermitian=True) @ np.ones(len(Sigma))
    return g / g.sum(), rank


# ---------------------------------------------------------------- pooled estimate
def pooled(psi, Phi, n_units=None, level=0.95, kappa=None):
    """EB pooled estimate with its valid (sandwich) and invalid (working) standard errors."""
    psi = np.asarray(psi, float)
    Sigma = Phi.T @ Phi
    v = np.diag(Sigma).copy()
    n_units = Phi.shape[0] if n_units is None else n_units
    out = {"J": len(psi)}
    t_pair, t_dl = tau2_pairwise(psi, v), tau2_dl(psi, v)
    t_thr, fired, ratio = tau2_thresholded(t_pair, v, n_units, kappa)
    out.update(tau_pair=np.sqrt(t_pair), tau_dl=np.sqrt(t_dl), tau2_over_mean_v=ratio,
               kappa=np.log(n_units) if kappa is None else kappa, threshold_fired=fired)
    out.update(dispersion_z=float(dispersion_z(psi, Sigma)))
    for name, t in (("pair", t_pair), ("thr", t_thr), ("dl", t_dl)):
        w = eb_weights(v, t)
        est = float(w @ psi)
        se = float(np.sqrt(w @ Sigma @ w))                 # observation-level sandwich
        se_work = float(np.sqrt(1 / np.sum(1 / (v + t))))  # working independence: NOT valid
        out.update({f"psi_eb_{name}": est, f"se_sandwich_{name}": se, f"se_working_{name}": se_work,
                    f"ci_lo_{name}": est - Z[level] * se, f"ci_hi_{name}": est + Z[level] * se})
    k = int(np.argmin(v))
    out.update(best_single_index=k, psi_best_single=float(psi[k]), se_best_single=float(np.sqrt(v[k])))
    wg, rank = gls_weights(Sigma)
    out.update(psi_gls=float(wg @ psi), se_gls=float(np.sqrt(max(wg @ Sigma @ wg, 0.0))),
               sigma_rank=rank, gls_min_weight=float(wg.min()))
    c = Sigma / np.sqrt(np.outer(v, v))
    iu = np.triu_indices(len(psi), 1)
    out.update(corr_mean=float(c[iu].mean()) if len(psi) > 1 else np.nan,
               corr_min=float(c[iu].min()) if len(psi) > 1 else np.nan,
               corr_max=float(c[iu].max()) if len(psi) > 1 else np.nan)
    return out


# ---------------------------------------------------------------- prediction intervals
def _conf_q(scores, p):
    m = len(scores)
    idx = int(np.floor(p * (m + 1))) if p <= 0.5 else int(np.ceil(p * (m + 1)))
    s = np.sort(scores)
    return np.nan if idx < 1 or idx > m else s[idx - 1]


def split_conformal(psi, v, v_new, alpha=None, B=2000, seed=42):
    """Split conformal interval of Section 6 (the object of Theorem 3), over B random splits.

    Training set of size floor(J/2), calibration set the rest. With alpha=None the finest
    two-sided level 1 - 2/(|C|+1) is used (interval between the extreme calibration scores).
    Returns the interval for the first split (a single analysis) and the median and spread
    of the endpoints over the B splits (stability of the split).
    """
    psi, v = np.asarray(psi, float), np.asarray(v, float)
    J = len(psi)
    n_tr = J // 2
    n_cal = J - n_tr
    alpha = 2 / (n_cal + 1) if alpha is None else alpha
    rng = np.random.default_rng(seed)
    los, his = [], []
    for _ in range(B):
        idx = rng.permutation(J)
        T, C = idx[:n_tr], idx[n_tr:]
        t = tau2_pairwise(psi[T], v[T])
        m = float(eb_weights(v[T], t) @ psi[T])
        S = (psi[C] - m) / np.sqrt(t + v[C])
        sc = np.sqrt(t + v_new)
        los.append(m + _conf_q(S, alpha / 2) * sc)
        his.append(m + _conf_q(S, 1 - alpha / 2) * sc)
    los, his = np.array(los), np.array(his)
    return dict(level=1 - alpha, n_train=n_tr, n_cal=n_cal, lo_first=los[0], hi_first=his[0],
                lo_median=float(np.nanmedian(los)), hi_median=float(np.nanmedian(his)),
                lo_sd=float(np.nanstd(los, ddof=1)), hi_sd=float(np.nanstd(his, ddof=1)), B=B)


def loo_conformal(psi, v, v_new):
    """Leave-one-out variant (uses all J scores; not covered by Theorem 3)."""
    psi, v = np.asarray(psi, float), np.asarray(v, float)
    J = len(psi)
    S = np.empty(J)
    for j in range(J):
        keep = np.arange(J) != j
        t = tau2_pairwise(psi[keep], v[keep])
        S[j] = (psi[j] - eb_weights(v[keep], t) @ psi[keep]) / np.sqrt(t + v[j])
    alpha = 2 / (J + 1)
    t = tau2_pairwise(psi, v)
    m = float(eb_weights(v, t) @ psi)
    sc = np.sqrt(t + v_new)
    return dict(level=1 - alpha, lo=m + _conf_q(S, alpha / 2) * sc, hi=m + _conf_q(S, 1 - alpha / 2) * sc,
                scores=S)
