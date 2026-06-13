## lalonde_functionals.R
## Computes J=15 HCT-framed identification functionals for the LaLonde data.
##
## Structure (Zhu et al. ECT framework):
##   theta_1 = E[Y(1) | S=1]  (NSW treated mean, fixed across all j)
##
##   RCT-only (j=1,2,3):
##     j=1: DIM      -- crude difference, NSW treated vs NSW control
##     j=2: IPW-NB   -- Hajek with e(X) = P(A=1|X, S=1), logistic on NSW
##     j=3: OR-NB    -- outcome regression fit on NSW controls, averaged over NSW
##
##   EC-augmented (j=4..15), for k in {CPS-2, CPS-3, PSID-2, PSID-3}:
##     j=3+3(k-1)+1: DIM-k  -- crude EC group mean
##     j=3+3(k-1)+2: OM-k   -- outcome model fit on NSW_ctrl + EC_k,
##                              predicted at NSW covariates
##     j=3+3(k-1)+3: IPW-k  -- sampling PS pi_k(X)=P(S=1|X,S in {1,k}),
##                              Hajek reweighting EC_k to NSW distribution
##
## Each psi_j = theta_1 - theta_0^(j),  v_j = Var(psi_j) via influence functions.
##
## EB estimator: psi_EB = sum_j w_j * psi_j / sum_j w_j,  w_j = 1/(v_j + tau2_hat)
## tau2_hat via marginal likelihood / pairwise method-of-moments

library(dplyr)

set.seed(42)

dir_data <- "~/Dropbox/Aplicaciones/Overleaf/shrinkage/data"
d <- read.csv(file.path(dir_data, "lalonde_combined.csv"))

covs <- c("age", "education", "black", "hispanic", "married", "nodegree", "re74", "re75")

# ── Partition ────────────────────────────────────────────────────────────────
nsw_trt  <- d |> filter(Z == 0, A == 1)   # NSW treated  (n=185)
nsw_ctrl <- d |> filter(Z == 0, A == 0)   # NSW control  (n=260)
nsw      <- d |> filter(Z == 0)           # full RCT

ec_list  <- list(
  "CPS-2"  = d |> filter(Z == 1),         # n=2369
  "CPS-3"  = d |> filter(Z == 2),         # n=429
  "PSID-2" = d |> filter(Z == 3),         # n=253
  "PSID-3" = d |> filter(Z == 4)          # n=128
)

n1 <- nrow(nsw_trt)
theta1 <- mean(nsw_trt$Y)

cat(sprintf("theta_1 (NSW treated mean re78): $%.0f  (n=%d)\n", theta1, n1))
cat(sprintf("theta_0 benchmark (NSW ctrl): $%.0f\n", mean(nsw_ctrl$Y)))

# ── Influence-function variance helpers ──────────────────────────────────────
var_mean <- function(y) var(y) / length(y)   # Var(sample mean)

# Hajek IPW variance: influence function phi_i = w_i * y_i / mean(w)
var_hajek <- function(y, w) {
  mu <- sum(w * y) / sum(w)
  phi <- w * (y - mu) / mean(w)
  var(phi) / length(phi)
}

# ── RCT-only functionals ─────────────────────────────────────────────────────

## j=1: DIM
theta0_1 <- mean(nsw_ctrl$Y)
v_0_1    <- var_mean(nsw_ctrl$Y)
psi_1    <- theta1 - theta0_1
v_1      <- var_mean(nsw_trt$Y) + v_0_1

## j=2: IPW-NB (trial propensity e(X) on NSW, Hajek)
fml_ps   <- as.formula(paste("A ~", paste(covs, collapse = " + ")))
fit_e    <- glm(fml_ps, data = nsw, family = binomial)
nsw_ctrl$e_hat <- predict(fit_e, newdata = nsw_ctrl, type = "response")
# Hajek: reweight NSW_ctrl by e/(1-e) to get E[Y(0)|S=1]
w2 <- nsw_ctrl$e_hat / (1 - nsw_ctrl$e_hat)
theta0_2 <- sum(w2 * nsw_ctrl$Y) / sum(w2)
v_0_2    <- var_hajek(nsw_ctrl$Y, w2)
psi_2    <- theta1 - theta0_2
v_2      <- var_mean(nsw_trt$Y) + v_0_2

## j=3: OR-NB (outcome model on NSW_ctrl, predict at NSW covariates)
fml_or  <- as.formula(paste("Y ~", paste(covs, collapse = " + ")))
fit_or  <- lm(fml_or, data = nsw_ctrl)
mu_hat3 <- predict(fit_or, newdata = nsw_trt)   # predict at treated covs
theta0_3 <- mean(mu_hat3)
# Variance of OR estimator via influence function residuals
resid3   <- nsw_ctrl$Y - predict(fit_or, newdata = nsw_ctrl)
v_0_3    <- var(resid3) / nrow(nsw_ctrl) + var(mu_hat3) / n1
psi_3    <- theta1 - theta0_3
v_3      <- var_mean(nsw_trt$Y) + v_0_3

# ── EC-augmented functionals ─────────────────────────────────────────────────

results <- data.frame(
  j      = integer(),
  label  = character(),
  psi    = numeric(),
  v      = numeric(),
  se     = numeric(),
  theta0 = numeric(),
  stringsAsFactors = FALSE
)

add_row <- function(j, lbl, psi, v, th0) {
  results <<- rbind(results, data.frame(
    j=j, label=lbl, psi=psi, v=v, se=sqrt(v), theta0=th0,
    stringsAsFactors=FALSE))
}

add_row(1, "DIM-RCT",    psi_1, v_1, theta0_1)
add_row(2, "IPW-NB",     psi_2, v_2, theta0_2)
add_row(3, "OR-NB",      psi_3, v_3, theta0_3)

ec_names <- names(ec_list)
for (ki in seq_along(ec_list)) {
  ec_k <- ec_list[[ki]]
  kname <- ec_names[ki]
  base_j <- 3 + 3 * (ki - 1)

  ## DIM-k: crude mean of EC group
  theta0_dim <- mean(ec_k$Y)
  v_0_dim    <- var_mean(ec_k$Y)
  psi_dim    <- theta1 - theta0_dim
  v_dim      <- var_mean(nsw_trt$Y) + v_0_dim
  add_row(base_j + 1, paste0("DIM-", kname), psi_dim, v_dim, theta0_dim)

  ## OM-k: outcome model on NSW_ctrl + EC_k, predict at NSW_trt covariates
  pool_om <- bind_rows(nsw_ctrl, ec_k)
  fit_om  <- lm(fml_or, data = pool_om)
  mu_om   <- predict(fit_om, newdata = nsw_trt)
  theta0_om <- mean(mu_om)
  resid_om  <- pool_om$Y - predict(fit_om, newdata = pool_om)
  v_0_om    <- var(resid_om) / nrow(pool_om) + var(mu_om) / n1
  psi_om    <- theta1 - theta0_om
  v_om      <- var_mean(nsw_trt$Y) + v_0_om
  add_row(base_j + 2, paste0("OM-", kname), psi_om, v_om, theta0_om)

  ## IPW-k: sampling PS pi_k(X) = P(S=1 | X, S in {NSW, EC_k})
  ##   Label: S=1 for NSW (both treated+ctrl), S=0 for EC_k
  nsw_s  <- nsw  |> mutate(S = 1L)
  ec_s   <- ec_k |> mutate(S = 0L)
  pool_s <- bind_rows(nsw_s, ec_s)
  fml_s  <- as.formula(paste("S ~", paste(covs, collapse = " + ")))
  fit_s  <- glm(fml_s, data = pool_s, family = binomial)
  ec_k$pi_hat <- predict(fit_s, newdata = ec_k, type = "response")
  # Hajek: w = pi/(1-pi) reweights EC_k to the NSW marginal
  w_ipw   <- ec_k$pi_hat / (1 - ec_k$pi_hat)
  theta0_ipw <- sum(w_ipw * ec_k$Y) / sum(w_ipw)
  v_0_ipw    <- var_hajek(ec_k$Y, w_ipw)
  psi_ipw    <- theta1 - theta0_ipw
  v_ipw      <- var_mean(nsw_trt$Y) + v_0_ipw
  add_row(base_j + 3, paste0("IPW-", kname), psi_ipw, v_ipw, theta0_ipw)
}

# ── Print raw results ─────────────────────────────────────────────────────────
cat("\n── Raw functionals ──────────────────────────────────────────────────────\n")
cat(sprintf("%-20s  %8s  %7s  %8s  %8s\n",
            "Functional", "psi_hat", "SE", "95% CI lo", "95% CI hi"))
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("%-20s  %8.0f  %7.0f  %8.0f  %8.0f\n",
              r$label, r$psi, r$se,
              r$psi - 1.96 * r$se,
              r$psi + 1.96 * r$se))
}

# ── Empirical Bayes shrinkage ─────────────────────────────────────────────────

## tau2 via marginal method-of-moments (Efron-Morris style):
##   E[sum (psi_j - psi_bar_w)^2 / v_j] ≈ (J-1) + tau2 * sum(1/v_j - sum(1/v_j^2)/sum(1/v_j))
J   <- nrow(results)
v   <- results$v
psi <- results$psi

# Precision-weighted grand mean
w_prec  <- 1 / v
psi_bar <- sum(w_prec * psi) / sum(w_prec)

# Method-of-moments for tau2
Q <- sum((psi - psi_bar)^2 / v)
C <- sum(w_prec) - sum(w_prec^2) / sum(w_prec)
tau2_hat_mm <- max(0, (Q - (J - 1)) / C)

cat(sprintf("\n── Empirical Bayes (method-of-moments) ──────────────────────────────────\n"))
cat(sprintf("Q statistic: %.2f  (df = J-1 = %d)\n", Q, J - 1))
cat(sprintf("tau2_hat (MoM): %.0f   [tau_hat = $%.0f]\n",
            tau2_hat_mm, sqrt(tau2_hat_mm)))

## EB posterior means
w_eb     <- 1 / (v + tau2_hat_mm)
psi_eb   <- sum(w_eb * psi) / sum(w_eb)
v_eb     <- 1 / sum(w_eb)

cat(sprintf("psi_EB: $%.0f   SE_EB: $%.0f\n", psi_eb, sqrt(v_eb)))
cat(sprintf("Experimental benchmark (DIM-RCT): $%.0f\n", results$psi[1]))

## Shrinkage factors
results$shrinkage <- (tau2_hat_mm) / (v + tau2_hat_mm)
results$psi_EB_j  <- psi_bar + results$shrinkage * (psi - psi_bar)

# ── Summary table ─────────────────────────────────────────────────────────────
cat("\n── Full summary (sorted by psi) ─────────────────────────────────────────\n")
cat(sprintf("%-20s  %8s  %7s  %8s  %6s\n",
            "Functional", "psi_hat", "SE", "psi_EB_j", "shrink%"))
ord <- order(results$psi)
for (i in ord) {
  r <- results[i, ]
  cat(sprintf("%-20s  %8.0f  %7.0f  %8.0f  %5.1f%%\n",
              r$label, r$psi, r$se, r$psi_EB_j,
              100 * r$shrinkage))
}

cat(sprintf("\n%-20s  %8.0f  %7.0f\n", "psi_EB (global)", psi_eb, sqrt(v_eb)))
cat(sprintf("%-20s  %8.0f\n", "Expt. benchmark", results$psi[1]))

# ── Sub-analysis: covariate-adjusted only (J=11, exclude crude DIM-k) ────────
cat("\n── Sub-analysis: J=11 covariate-adjusted functionals only ──────────────\n")
cat("   (RCT-only 3 + OM×4 + IPW×4; excludes crude DIM-k with outcome drift)\n\n")

adj_idx <- !grepl("^DIM-CPS|^DIM-PSID", results$label)
r11  <- results[adj_idx, ]
psi11 <- r11$psi
v11   <- r11$v
J11   <- nrow(r11)
w11   <- 1 / v11
pbar11 <- sum(w11 * psi11) / sum(w11)
Q11   <- sum((psi11 - pbar11)^2 / v11)
C11   <- sum(w11) - sum(w11^2) / sum(w11)
tau2_11 <- max(0, (Q11 - (J11 - 1)) / C11)
web11  <- 1 / (v11 + tau2_11)
peb11  <- sum(web11 * psi11) / sum(web11)
veb11  <- 1 / sum(web11)

cat(sprintf("Q statistic: %.2f  (df = J-1 = %d)\n", Q11, J11 - 1))
cat(sprintf("tau2_hat (MoM): %.0f   [tau_hat = $%.0f]\n", tau2_11, sqrt(tau2_11)))
cat(sprintf("psi_EB (J=11): $%.0f   SE_EB: $%.0f\n", peb11, sqrt(veb11)))
cat(sprintf("Expt. benchmark: $%.0f\n", results$psi[1]))

## RCT-only sub-analysis (J=3)
cat("\n── Sub-analysis: RCT-only J=3 (tau2=0 regime) ───────────────────────────\n")
r3  <- results[1:3, ]
psi3 <- r3$psi; v3 <- r3$v
w3   <- 1/v3; pbar3 <- sum(w3*psi3)/sum(w3)
Q3   <- sum((psi3-pbar3)^2/v3)
tau2_3 <- max(0, (Q3-2)/(sum(w3)-sum(w3^2)/sum(w3)))
peb3 <- pbar3   # tau2=0 so EB = precision-weighted mean
veb3 <- 1/sum(w3)
cat(sprintf("Q statistic: %.2f  (df = 2)\n", Q3))
cat(sprintf("tau2_hat (MoM): %.0f  (%.0f  if > 0)\n", tau2_3, tau2_3))
cat(sprintf("psi_EB (J=3, tau2=0): $%.0f   SE_EB: $%.0f\n", peb3, sqrt(veb3)))

results$adj11_psi_EB <- NA_real_
results$adj11_psi_EB[adj_idx] <- peb11 + (tau2_11/(v11+tau2_11)) * (psi11 - peb11)

# ── Save ─────────────────────────────────────────────────────────────────────
out <- results |>
  select(j, label, psi, v, se, theta0, shrinkage, psi_EB_j) |>
  mutate(tau2_hat = tau2_hat_mm, psi_EB = psi_eb, v_EB = v_eb,
         theta1 = theta1)

write.csv(out, file.path(dir_data, "lalonde_functionals.csv"), row.names = FALSE)
cat(sprintf("\nSaved: data/lalonde_functionals.csv  (%d functionals)\n", J))
