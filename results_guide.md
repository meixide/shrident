# Reading the HAPC LaLonde results

All five output CSVs are in `data/` (working copy) and permanently in
`$STORE/shrinkage_results/` (backed up). They were produced on 2026-06-14 by a
full-data run (`HAPC_FRAC=1.0`) on CESGA FT3. All 10 CV fits returned
**INTERIOR**, so the regularisation grid was adequate.

---

## `lalonde_functionals_hapc.csv` — the 15 identifying functionals

One row per functional (j=1..15). Columns:

| Column | Meaning |
|---|---|
| `j` | Functional index (matches the theory paper's numbering) |
| `label` | `{estimator}-{dataset}`: estimator ∈ {DIM, OM, IPW, OR}, dataset ∈ {RCT, NB, CPS-2, CPS-3, PSID-2, PSID-3} |
| `psi` | Point estimate (ATT in $, same units as LaLonde 1986) |
| `v` | Variance estimate (not SE²; see note below) |
| `se` | Standard error = √v |
| `theta0` | Plug-in null (used as the shrinkage target; computed from the same fit) |

**Notes on `v`:** this is the semiparametric variance of the functional, not a
sandwich. For IPW functionals it can be large (IPW-CPS-2 SE ≈ 1618) because
the propensity scores are non-trivial to stabilise.

**Key numbers to report:**
- RCT benchmark (j=1, DIM-RCT): **ψ ≈ 1794** (se 671) — the ground truth.
- DIM on non-experimental comparisons (j=4, j=7, j=10, j=13) are badly biased
  (−3822, −635, −3647, +1070) — these demonstrate the selection problem.
- OM and IPW on non-experimental data recover estimates in the **1384–1839**
  range, consistent with the RCT benchmark.

---

## `lalonde_eb_hapc.csv` — empirical Bayes shrinkage

Three pooling choices; each has a moment-matching (`_mom`) and a pairwise
(`_pair`) shrinkage variant.

| Column | Meaning |
|---|---|
| `analysis` | Which subset of the J functionals was pooled |
| `J` | Number of functionals in the pool |
| `tau_mom` / `tau_pair` | Estimated heterogeneity τ (0 = homogeneous pool) |
| `psi_eb_mom` / `psi_eb_pair` | EB pooled ATT estimate |
| `se_eb_mom` / `se_eb_pair` | SE of the EB estimate |

**Key numbers:**
- RCT-only pool (J=3): **ψ_EB ≈ 1770** (se 389) — tight, no heterogeneity (τ=0).
- Covariate-adjusted pool (J=11, excludes the three biased DIMs): **ψ_EB ≈ 1552** (se 228) — the headline number for the paper; narrowest SE.
- Full pool (J=15, includes the biased DIMs): τ > 0, estimates around **590–656** — shrinkage towards the biased functionals pulls the estimate down; illustrates the cost of pooling with misspecified comparisons.

---

## `lalonde_loo_hapc.csv` — leave-one-out predictive scores

Same 15 rows as `lalonde_functionals_hapc.csv`, plus:

| Column | Meaning |
|---|---|
| `S_loo` | LOO predictive score for functional j (higher = better supported by the other 14) |

Positive S_loo means the held-out functional is well-predicted by the pool.
Negative values (j=4 DIM-CPS-2: −3.63; j=10 DIM-PSID-2: −2.94; j=7 DIM-CPS-3:
−0.72) flag the biased DIM functionals as outliers — they are inconsistent with
the rest of the evidence. This is the formal diagnostic for "this functional
should not be in the pool."

---

## `lalonde_pi_hapc.csv` — prediction intervals for a new functional

Three interval methods for the marginal predictive distribution of a new
identifying functional drawn from the same population:

| Column | Meaning |
|---|---|
| `method` | Interval type and nominal coverage |
| `lower` / `upper` | Endpoints in $ |
| `width` | `upper − lower` |

The Normal 90% interval (−2497, 3671) is the parametric baseline. The
split-conformal (77.8%) and LOO-conformal (87.5%) intervals are
distribution-free; their nominal coverages differ from 90% because of the small
J. All three bracket roughly the same range.

---

## `lalonde_functionals_ols_sub.csv` — OLS/logistic baseline

Same schema as `lalonde_functionals_hapc.csv`. The DIM rows are identical (no
model involved); the OM and IPW rows use OLS (gaussian) or logistic regression
instead of HAPC. Use this for the apples-to-apples comparison table in the
paper: HAPC vs. linear base learners on the same full sample.

Notable differences vs. HAPC:
- IPW-CPS-2 OLS: ψ ≈ 1014 (se 688) vs. HAPC 1451 (se 1618) — OLS propensity is
  better stabilised here but the point estimate is lower.
- IPW-PSID-3 OLS: ψ ≈ 2001 (se 1411) vs. HAPC 1220 (se 1209) — OLS overshoots.
- OM estimates are close across methods (within ~100 $), consistent with the
  outcome surfaces being smooth.

---

## Quick Python snippet to load everything

```python
import pandas as pd, pathlib

D = pathlib.Path("data")   # or pathlib.Path(os.environ["STORE"]) / "shrinkage_results"

psi   = pd.read_csv(D / "lalonde_functionals_hapc.csv").set_index("j")
eb    = pd.read_csv(D / "lalonde_eb_hapc.csv")
loo   = pd.read_csv(D / "lalonde_loo_hapc.csv").set_index("j")
pi    = pd.read_csv(D / "lalonde_pi_hapc.csv")
ols   = pd.read_csv(D / "lalonde_functionals_ols_sub.csv").set_index("j")

# headline EB estimate (covariate-adjusted pool, pairwise shrinkage)
row = eb[eb.analysis.str.startswith("Covariate")].iloc[0]
print(f"EB ATT = {row.psi_eb_pair:.0f}  (SE {row.se_eb_pair:.0f})")
```
