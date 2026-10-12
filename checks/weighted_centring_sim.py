"""Consistency needs the weighted discrepancy to vanish: which weights, under which pattern.

psi_hat_j = psi* + eps_j + xi_j, psi* = 0, j = 1..J. Sampling errors xi ~ N(0, Omega/n) with
Omega = rho sigma sigma' + (1-rho) diag(sigma^2), rho = 0.2, sigma_j iid U[0.6, 2], n = 400.
Each formula also has a design weight p_j (Gamma(2) draws, normalised), e.g. a stratum share.
Discrepancy patterns, all with standard deviation about tau = 0.1:
  balanced           eps_j iid N(0, tau^2), independent of sigma_j and p_j
  tied to precision  eps_j = tau * standardised sigma_j: mean zero, predictable from precision
  common sign        eps_j iid N(tau, tau^2)
  design-centred     sum_j p_j eps_j = 0 exactly, eps_j increasing in p_j (strata of one population)
Weights: uniform; precision w_j(tau2-hat) with the pairwise estimator; design p_j.
Reported: bias and root mean squared error of the pooled estimate over 4000 replications.
The sampling error alone gives an RMSE of about sqrt(rho) E(sigma) / sqrt(n) = 0.03.

Run from the repo root:  python checks/weighted_centring_sim.py
"""
import numpy as np

TAU, RHO, N = 0.1, 0.2, 400


def pairwise_tau2(psi, v):
    return max(0.0, psi.var(ddof=1) - v.mean())      # identical to the average over pairs


def draw(pattern, J, rng):
    sig = rng.uniform(0.6, 2.0, J)
    p = rng.gamma(2.0, size=J)
    p /= p.sum()
    if pattern == "balanced":
        eps = rng.normal(0, TAU, J)
    elif pattern == "tied to precision":
        eps = TAU * (sig - 1.3) / (1.4 / np.sqrt(12))
    elif pattern == "common sign":
        eps = rng.normal(TAU, TAU, J)
    else:                                             # design-centred
        e = TAU * (J * p - 1) / np.sqrt(0.5) + rng.normal(0, TAU / 2, J)
        eps = e - p @ e
    f = rng.standard_normal()
    xi = sig * (np.sqrt(RHO) * f + np.sqrt(1 - RHO) * rng.standard_normal(J)) / np.sqrt(N)
    return eps + xi, sig ** 2 / N, p


def main(reps=4000, seed=11):
    rng = np.random.default_rng(seed)
    print(f"{'pattern':18s} {'J':>4s} | " + " | ".join(f"{w:>17s}" for w in ("uniform", "precision", "design")))
    print(f"{'':18s} {'':>4s} | " + " | ".join(f"{'bias':>8s} {'rmse':>8s}" for _ in range(3)))
    for pattern in ("balanced", "tied to precision", "common sign", "design-centred"):
        for J in (10, 40, 160):
            out = {k: [] for k in ("uniform", "precision", "design")}
            for _ in range(reps):
                psi, v, p = draw(pattern, J, rng)
                wp = 1 / (v + pairwise_tau2(psi, v))
                out["uniform"].append(psi.mean())
                out["precision"].append(wp @ psi / wp.sum())
                out["design"].append(p @ psi)
            cells = [f"{np.mean(x):8.3f} {np.sqrt(np.mean(np.square(x))):8.3f}" for x in out.values()]
            print(f"{pattern:18s} {J:4d} | " + " | ".join(cells))


if __name__ == "__main__":
    main()
