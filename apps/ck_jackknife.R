## ck_jackknife.R
## Card-Krueger application: restaurant-level influence contributions of the 25
## estimators, by a delete-one jackknife stratified by state.
##
##   Rscript apps/ck_jackknife.R results/ck94_jackknife
##
## The estimators are those of data/ck94_functionals.R (difference in means, outcome
## model, weighting; nine Pennsylvania subgroups); the script stops with status 3 if
## its full-sample estimates differ from data/ck94_functionals.csv. It then deletes one
## restaurant at a time, recomputes the 25 estimates, and writes
##   psi.csv   label, estimate, jackknife variance and the published variance
##   Phi.csv   one row per restaurant, one column per estimator, scaled so that
##             Sigma = Phi'Phi is the stratified jackknife covariance matrix:
##             Phi[i, j] = sqrt((n_s - 1) / n_s) * (mean_s psi_j(-.) - psi_j(-i)),
##             s the state of restaurant i
##   MANIFEST
## An adjusted estimator whose covariates do not vary in a leave-one-out sample equals
## the difference in means there; the number of such cases is recorded.

suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript apps/ck_jackknife.R <output directory>")
out <- args[1]
if (dir.exists(out)) stop(out, " exists. Choose a new name; a run directory is never reused.")
dir.create(out, recursive = TRUE)

raw <- read.csv("data/ck94.csv") |>
  filter(!is.na(delta_fte)) |>
  mutate(treat = as.integer(state == 1),
         chain_bk = as.integer(chain == 1), chain_kfc = as.integer(chain == 2),
         chain_rr = as.integer(chain == 3), unit = row_number())
cov_cols <- c("chain_bk", "chain_kfc", "chain_rr", "co_owned")
groups <- list("PA-all" = function(p) p, "PA-north" = function(p) filter(p, pa1 == 1),
               "PA-east" = function(p) filter(p, pa2 == 1), "PA-BK" = function(p) filter(p, chain == 1),
               "PA-KFC" = function(p) filter(p, chain == 2), "PA-RR" = function(p) filter(p, chain == 3),
               "PA-Wendys" = function(p) filter(p, chain == 4), "PA-corp" = function(p) filter(p, co_owned == 1),
               "PA-fran" = function(p) filter(p, co_owned == 0))
safe_cols <- function(X) X[, apply(X, 2, function(x) var(x) > 0), drop = FALSE]

# The three estimators, as in data/ck94_functionals.R (estimates only). NA = not defined.
psi_dim <- function(nj, pa_k) mean(nj$delta_fte) - mean(pa_k$delta_fte)
psi_om <- function(nj, pa_k) {
  X_pa <- safe_cols(pa_k[, cov_cols, drop = FALSE])
  if (ncol(X_pa) == 0 || nrow(pa_k) <= ncol(X_pa) + 3) return(NA_real_)
  fit <- lm(y ~ ., data = data.frame(y = pa_k$delta_fte, X_pa))
  mean(nj$delta_fte) - mean(suppressWarnings(predict(fit, newdata = nj[, colnames(X_pa), drop = FALSE])))
}
psi_ipw <- function(nj, pa_k) {
  if (nrow(pa_k) < 10) return(NA_real_)
  X_pa <- safe_cols(pa_k[, cov_cols, drop = FALSE])
  if (ncol(X_pa) == 0) return(NA_real_)
  X_nj <- nj[, colnames(X_pa), drop = FALSE]
  df <- data.frame(A = c(rep(1L, nrow(X_nj)), rep(0L, nrow(X_pa))), rbind(X_nj, X_pa))
  fit <- suppressWarnings(glm(A ~ ., data = df, family = binomial))
  ps <- pmax(pmin(suppressWarnings(predict(fit, newdata = data.frame(X_pa), type = "response")), 0.99), 0.01)
  w <- ps / (1 - ps)
  mean(nj$delta_fte) - sum(w * pa_k$delta_fte) / sum(w)
}
all_psi <- function(d) {
  nj <- filter(d, treat == 1); pa <- filter(d, treat == 0)
  r <- c()
  for (nm in names(groups)) {
    pa_k <- groups[[nm]](pa)
    r[paste0("DIM-", nm)] <- psi_dim(nj, pa_k)
    r[paste0("OM-", nm)]  <- psi_om(nj, pa_k)
    r[paste0("IPW-", nm)] <- psi_ipw(nj, pa_k)
  }
  r
}

# ---- full sample, checked against the published table
published <- read.csv("data/ck94_functionals.csv")
full <- all_psi(raw)
defined <- names(full)[!is.na(full)]
if (!setequal(defined, published$label)) stop("the set of estimators differs from data/ck94_functionals.csv")
psi <- full[published$label]
gap <- max(abs(psi - published$psi))
cat(sprintf("full-sample estimates: %d, largest gap to data/ck94_functionals.csv %.3g\n", length(psi), gap))
if (gap > 1e-8) { cat("CHECK FAILED\n"); quit(status = 3) }

# ---- delete-one jackknife, stratified by state
n <- nrow(raw)
loo <- matrix(NA_real_, n, length(psi), dimnames = list(NULL, names(psi)))
fallback <- setNames(integer(length(psi)), names(psi))
for (i in seq_len(n)) {
  r <- all_psi(raw[-i, ])[names(psi)]
  miss <- names(r)[is.na(r)]
  for (m in miss) {            # no covariate varies: the adjusted estimator is the difference in means
    r[m] <- r[sub("^(OM|IPW)-", "DIM-", m)]
    fallback[m] <- fallback[m] + 1L
  }
  loo[i, ] <- r
}
Phi <- loo
for (s in c(0, 1)) {
  k <- raw$treat == s
  ns <- sum(k)
  Phi[k, ] <- sqrt((ns - 1) / ns) * sweep(-loo[k, , drop = FALSE], 2, -colMeans(loo[k, , drop = FALSE]))
}
Sigma <- crossprod(Phi)
res <- data.frame(label = names(psi), ctrl = published$ctrl, method = published$method, n_ctrl = published$n_ctrl,
                  psi = as.numeric(psi), v = diag(Sigma), se = sqrt(diag(Sigma)),
                  v_published = published$v, se_published = published$se, loo_fallbacks = as.integer(fallback))
write.csv(res, file.path(out, "psi.csv"), row.names = FALSE)
write.csv(data.frame(unit = raw$unit, state = ifelse(raw$treat == 1, "NJ", "PA"), Phi, check.names = FALSE),
          file.path(out, "Phi.csv"), row.names = FALSE)
writeLines(c(
  paste0("run=", basename(out)),
  paste0("created_utc=", format(Sys.time(), "%Y-%m-%dT%H:%M:%SZ", tz = "UTC")),
  paste0("git_commit=", tryCatch(system("git rev-parse HEAD", intern = TRUE), error = function(e) "unavailable")),
  paste0("input=data/ck94.csv md5=", tools::md5sum("data/ck94.csv")),
  paste0("restaurants=", n, " NJ=", sum(raw$treat == 1), " PA=", sum(raw$treat == 0)),
  paste0("estimators=", length(psi)),
  "variance=delete-one jackknife, stratified by state, factor (n_s-1)/n_s",
  paste0("full_sample_max_gap_to_published=", signif(gap, 3)),
  paste0("loo_fallbacks_total=", sum(fallback), " (adjusted estimator set to the difference in means when no covariate varies)"),
  paste0("R_version=", R.version.string)), file.path(out, "MANIFEST"))
print(res[, c("label", "psi", "se", "se_published", "loo_fallbacks")], digits = 4)
