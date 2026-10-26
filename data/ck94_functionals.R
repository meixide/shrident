## ck94_functionals.R
## Card & Krueger (1994) DiD with multiple PA control groups.
## Computes J identification functionals for the ATT (NJ vs PA subsets),
## applies EB shrinkage and LOO conformal scores.
##
## Output: data/ck94_functionals.csv

library(dplyr)
library(sandwich)

dir_data <- "data"  # relative to the repository root

# ── Load and clean ─────────────────────────────────────────────────────────────
raw <- read.csv(file.path(dir_data, "ck94.csv")) |>
  filter(!is.na(delta_fte)) |>
  mutate(
    treat     = as.integer(state == 1),
    chain_bk  = as.integer(chain == 1),
    chain_kfc = as.integer(chain == 2),
    chain_rr  = as.integer(chain == 3)
    # chain_w (Wendy's, chain==4) is the reference category in regressions
  )

# NJ (treated) and PA (potential controls)
nj <- raw |> filter(treat == 1)
pa <- raw |> filter(treat == 0)

# Covariates available in BOTH NJ and PA (chain type and ownership)
cov_cols <- c("chain_bk", "chain_kfc", "chain_rr", "co_owned")

# ── Variance helpers ───────────────────────────────────────────────────────────

# HC2 variance for difference of means
var_dim <- function(y1, y0) var(y1) / length(y1) + var(y0) / length(y0)

# OM influence-function variance:
#   v = S²_NJ/n_NJ  +  (1/n_NJ) * sum over NJ of (ŷ_{0j} - ȳ̂_0)²
# (sampling variance of the predicted counterfactual mean)
var_om <- function(y_nj, yhat0_nj, r_pa, X_pa, fit) {
  n1 <- length(y_nj)
  # HC1 sandwich on PA regression for plug-in robustness
  vcov_pa <- tryCatch(vcovHC(fit, type = "HC1"), error = function(e) vcov(fit))
  v_nj    <- var(y_nj) / n1
  # Variance from prediction uncertainty: Var(x_i' β̂) averaged over NJ
  Xnj_mat <- model.matrix(~ chain_bk + chain_kfc + chain_rr + co_owned,
                           data = data.frame(X_pa)[rep(1, n1), , drop = FALSE])
  # Use NJ design matrix if available
  v_pred  <- mean(diag(crossprod(t(Xnj_mat), vcov_pa) * t(Xnj_mat)))
  v_nj + v_pred
}

# Hajek IPW variance for the PA mean component:
#   V̂_Hajek = (1/n0) * Σ (ω_i * (y_i - μ̂_w))²  where ω_i = w_i / Σw_j
var_ipw_pa <- function(y_pa, w_pa) {
  w_n   <- w_pa / sum(w_pa)
  mu_w  <- sum(w_n * y_pa)
  mean((w_n * (y_pa - mu_w))^2)
}

# ── Estimators ─────────────────────────────────────────────────────────────────

# Drop columns with zero variance (avoid rank deficiency in small subgroups)
safe_cols <- function(X) X[, apply(X, 2, function(x) var(x) > 0), drop = FALSE]

# Subset of NJ covariates to match a PA covariate subset
nj_covs <- function(pa_X) nj[, colnames(pa_X), drop = FALSE]

compute_dim <- function(pa_k) {
  y_nj <- nj$delta_fte; y_pa <- pa_k$delta_fte
  psi  <- mean(y_nj) - mean(y_pa)
  v    <- var_dim(y_nj, y_pa)
  data.frame(psi = psi, v = v, se = sqrt(v))
}

compute_om <- function(pa_k) {
  y_nj <- nj$delta_fte
  y_pa <- pa_k$delta_fte
  X_pa <- safe_cols(pa_k[, cov_cols, drop = FALSE])
  if (ncol(X_pa) == 0 || nrow(pa_k) <= ncol(X_pa) + 3) return(NULL)
  X_nj <- nj_covs(X_pa)
  df_pa <- data.frame(y = y_pa, X_pa)
  df_nj <- data.frame(X_nj);  colnames(df_nj) <- colnames(X_pa)
  fit   <- tryCatch(lm(y ~ ., data = df_pa), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  yhat0 <- predict(fit, newdata = df_nj)
  psi   <- mean(y_nj) - mean(yhat0)
  # Influence-function SE
  n1    <- length(y_nj)
  n0    <- length(y_pa)
  X_nj_mat <- model.matrix(~., data = df_nj)
  vcov_pa  <- tryCatch(vcovHC(fit, type = "HC1"), error = function(e) vcov(fit))
  # Var(ȳ_NJ) + Var(x̄_NJ' β̂)
  v_mean_nj <- var(y_nj) / n1
  xbar      <- colMeans(X_nj_mat)
  v_pred    <- as.numeric(t(xbar) %*% vcov_pa %*% xbar) / n1
  v <- v_mean_nj + v_pred
  if (!is.finite(v) || v <= 0) return(NULL)
  data.frame(psi = psi, v = v, se = sqrt(v))
}

compute_ipw <- function(pa_k) {
  if (nrow(pa_k) < 10) return(NULL)
  y_nj <- nj$delta_fte
  y_pa <- pa_k$delta_fte
  X_pa <- safe_cols(pa_k[, cov_cols, drop = FALSE])
  if (ncol(X_pa) == 0) return(NULL)
  X_nj <- nj_covs(X_pa)
  n1   <- nrow(X_nj); n0 <- nrow(X_pa)
  X_pool <- rbind(X_nj, X_pa)
  A_pool <- c(rep(1L, n1), rep(0L, n0))
  df_ps  <- data.frame(A = A_pool, X_pool)
  fit_ps <- tryCatch(glm(A ~ ., data = df_ps, family = binomial),
                     error = function(e) NULL)
  if (is.null(fit_ps)) return(NULL)
  ps_pa  <- predict(fit_ps, newdata = data.frame(X_pa), type = "response")
  # Guard against extreme propensities
  ps_pa  <- pmax(pmin(ps_pa, 0.99), 0.01)
  w_pa   <- ps_pa / (1 - ps_pa)
  mu_pa  <- sum(w_pa * y_pa) / sum(w_pa)
  psi    <- mean(y_nj) - mu_pa
  # Variance: Var(Ȳ_NJ) + Hajek variance of weighted PA mean
  w_n    <- w_pa / sum(w_pa)
  v_pa   <- mean((w_n * (y_pa - mu_pa))^2)   # Hajek
  v      <- var(y_nj) / n1 + v_pa
  if (!is.finite(v) || v <= 0) return(NULL)
  data.frame(psi = psi, v = v, se = sqrt(v))
}

# ── Define control group subsets ───────────────────────────────────────────────
ctrl_groups <- list(
  "PA-all"    = pa,
  "PA-north"  = pa |> filter(pa1 == 1),
  "PA-east"   = pa |> filter(pa2 == 1),
  "PA-BK"     = pa |> filter(chain == 1),
  "PA-KFC"    = pa |> filter(chain == 2),
  "PA-RR"     = pa |> filter(chain == 3),
  "PA-Wendys" = pa |> filter(chain == 4),
  "PA-corp"   = pa |> filter(co_owned == 1),
  "PA-fran"   = pa |> filter(co_owned == 0)
)

# ── Compute all functionals ───────────────────────────────────────────────────
rows <- list()
for (nm in names(ctrl_groups)) {
  pa_k <- ctrl_groups[[nm]]
  n_k  <- nrow(pa_k)

  # DIM-k
  r <- compute_dim(pa_k)
  rows[[paste0("DIM-", nm)]] <- data.frame(
    label = paste0("DIM-", nm), ctrl = nm, method = "DIM",
    n_ctrl = n_k, r)

  # OM-k
  r <- compute_om(pa_k)
  if (!is.null(r)) rows[[paste0("OM-", nm)]] <- data.frame(
    label = paste0("OM-", nm), ctrl = nm, method = "OM",
    n_ctrl = n_k, r)

  # IPW-k
  r <- compute_ipw(pa_k)
  if (!is.null(r)) rows[[paste0("IPW-", nm)]] <- data.frame(
    label = paste0("IPW-", nm), ctrl = nm, method = "IPW",
    n_ctrl = n_k, r)
}

d_res <- bind_rows(rows)
cat(sprintf("Computed %d functionals\n", nrow(d_res)))
print(d_res[, c("label", "psi", "se", "n_ctrl")])

# ── EB helpers ─────────────────────────────────────────────────────────────────
mom_tau2 <- function(psi, v) {
  w    <- 1 / v
  pbar <- sum(w * psi) / sum(w)
  Q    <- sum(w * (psi - pbar)^2)
  C    <- sum(w) - sum(w^2) / sum(w)
  max(0, (Q - (length(psi) - 1)) / C)
}

pairwise_tau2 <- function(psi, v) {
  pairs <- combn(length(psi), 2)
  vals  <- ((psi[pairs[1,]] - psi[pairs[2,]])^2 -
              (v[pairs[1,]] + v[pairs[2,]])) / 2
  max(0, mean(vals))
}

eb_posterior <- function(psi, v, tau2) {
  web  <- 1 / (v + tau2)
  psiEB <- sum(web * psi) / sum(web)
  list(psiEB = psiEB,
       psi_EB_j = psiEB + tau2 / (v + tau2) * (psi - psiEB))
}

cochran_Q <- function(psi, v) {
  w    <- 1 / v
  pbar <- sum(w * psi) / sum(w)
  sum(w * (psi - pbar)^2)
}

# ── All-J EB ──────────────────────────────────────────────────────────────────
J   <- nrow(d_res)
psi <- d_res$psi
v   <- d_res$v
tau2_all <- mom_tau2(psi, v)
eb_all   <- eb_posterior(psi, v, tau2_all)
Q_all    <- cochran_Q(psi, v)

d_res <- d_res |>
  mutate(
    tau2_hat = tau2_all,
    psiEB    = eb_all$psiEB,
    psi_EB_j = eb_all$psi_EB_j,
    lo95     = psi - 1.96 * se,
    hi95     = psi + 1.96 * se
  )

cat(sprintf("\nAll J=%d: Q=%.2f  crit(df=%d)=%.2f  tau=%.3f  psi_EB=%.3f\n",
    J, Q_all, J-1, qchisq(0.95, J-1), sqrt(tau2_all), eb_all$psiEB))

# ── Sub-analyses ──────────────────────────────────────────────────────────────
for (meth in c("DIM", "OM", "IPW")) {
  idx <- d_res$method == meth
  if (sum(idx) < 2) next
  psi_m <- psi[idx]; v_m <- v[idx]
  J_m   <- sum(idx)
  tau2_m <- mom_tau2(psi_m, v_m)
  eb_m   <- eb_posterior(psi_m, v_m, tau2_m)
  Q_m    <- cochran_Q(psi_m, v_m)
  cat(sprintf("  %-3s J=%d: Q=%.2f  crit(df=%d)=%.2f  tau=%.3f  psi_EB=%.3f\n",
      meth, J_m, Q_m, J_m-1, qchisq(0.95, J_m-1), sqrt(tau2_m), eb_m$psiEB))
}

# ── LOO conformal scores ───────────────────────────────────────────────────────
conformal_quantile <- function(scores, p) {
  Jc  <- length(scores)
  idx <- if (p <= 0.5) floor(p * (Jc + 1)) else ceiling(p * (Jc + 1))
  s   <- sort(scores)
  if (idx < 1 || idx > Jc) return(NA_real_)
  s[idx]
}

S_loo <- numeric(J)
tau2_loo <- numeric(J)
for (j in seq_len(J)) {
  tr  <- setdiff(seq_len(J), j)
  t2  <- pairwise_tau2(psi[tr], v[tr])
  eb  <- eb_posterior(psi[tr], v[tr], t2)
  tau2_loo[j] <- t2
  S_loo[j]    <- (psi[j] - eb$psiEB) / sqrt(t2 + v[j])
}

alpha_loo <- 2 / (J + 1)
q_lo <- conformal_quantile(S_loo, alpha_loo / 2)
q_hi <- conformal_quantile(S_loo, 1 - alpha_loo / 2)

tau2_pw  <- pairwise_tau2(psi, v)
v_new    <- median(v)
pi_lo    <- eb_all$psiEB + q_lo * sqrt(tau2_pw + v_new)
pi_hi    <- eb_all$psiEB + q_hi * sqrt(tau2_pw + v_new)
cat(sprintf("\nLOO conformal (J=%d, 1-alpha=%.1f%%): q_lo=%.3f  q_hi=%.3f\n",
    J, 100*(1-alpha_loo), q_lo, q_hi))
cat(sprintf("  PI (v_new=median): [%.3f, %.3f] FTE\n", pi_lo, pi_hi))

d_res$S_loo     <- S_loo
d_res$tau2_loo  <- tau2_loo

write.csv(d_res, file.path(dir_data, "ck94_functionals.csv"), row.names = FALSE)
cat("\nSaved ck94_functionals.csv\n")
