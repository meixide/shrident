## Cross-check of the Card-Krueger jackknife standard errors (apps/ck_jackknife.R):
## (1) a bootstrap stratified by state, 2000 draws, estimators as in the jackknife script;
## (2) the heteroskedasticity-consistent formula for the outcome-model estimator against
##     all of Pennsylvania.
## Run from the repo root:  Rscript checks/ck_jackknife_crosscheck.R > checks/ck_jackknife_crosscheck.out
suppressMessages({library(dplyr); library(sandwich)})
src <- readLines("apps/ck_jackknife.R")
a <- grep("^raw <- read.csv", src); b <- grep("^# ---- full sample", src)
eval(parse(text = src[a:(b - 1)]))            # data and estimator definitions only
full <- all_psi(raw); lab <- names(full)[!is.na(full)]
set.seed(1); B <- 2000
nj_i <- which(raw$treat == 1); pa_i <- which(raw$treat == 0)
bt <- matrix(NA_real_, B, length(lab), dimnames = list(NULL, lab))
for (r in 1:B) {
  d <- raw[c(sample(nj_i, replace = TRUE), sample(pa_i, replace = TRUE)), ]
  bt[r, ] <- tryCatch(all_psi(d)[lab], error = function(e) rep(NA_real_, length(lab)))
}
jk <- read.csv("results/ck94_jackknife_stratified/psi.csv")
out <- data.frame(label = lab, jackknife = jk$se, bootstrap_sd = apply(bt, 2, sd, na.rm = TRUE),
                  published = jk$se_published, draws_undefined = colSums(is.na(bt)))
print(out, digits = 3, row.names = FALSE)
nj <- filter(raw, treat == 1); pa <- filter(raw, treat == 0)
fit <- lm(delta_fte ~ chain_bk + chain_kfc + chain_rr + co_owned, data = pa)
X <- model.matrix(~ chain_bk + chain_kfc + chain_rr + co_owned, data = nj)
cat(sprintf("\nOM-PA-all, heteroskedasticity-consistent formula: %.3f\n",
            sqrt(var(nj$delta_fte - X %*% coef(fit)) / nrow(nj) + t(colMeans(X)) %*% vcovHC(fit, type = "HC1") %*% colMeans(X))))
cat(sprintf("ratio published / jackknife standard error, adjusted estimators: %.2f to %.2f\n",
            min((jk$se_published / jk$se)[jk$method != "DIM"]), max((jk$se_published / jk$se)[jk$method != "DIM"])))
