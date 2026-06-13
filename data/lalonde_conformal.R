## lalonde_conformal.R
## Split conformal prediction intervals for the J=15 LaLonde functionals,
## exactly following Section 5 of the paper (Theorem thm:conformal).
##
## Targets: psi_{J+1}, the latent target of a NEW identification functional
##   drawn from the same population. NOT a CI for psi*.
##
## Three analyses:
##   A. J=15 (tau2 > 0) -- conformal valid asymptotically
##   B. J=11 covariate-adjusted (tau2 = 0) -- conformal invalid; for contrast
##   C. LOO diagnostic: jackknife scores to check calibration

library(dplyr)

set.seed(42)
dir_data  <- "~/Dropbox/Aplicaciones/Overleaf/shrinkage/data"
dir_plots <- "~/Dropbox/Aplicaciones/Overleaf/shrinkage/plots"

d <- read.csv(file.path(dir_data, "lalonde_functionals.csv"))

# ── Helpers ───────────────────────────────────────────────────────────────────

## Pairwise tau2 estimator (Section 5, training step)
tau2_pairwise <- function(psi, v) {
  J  <- length(psi)
  pairs <- combn(J, 2)
  vals  <- ((psi[pairs[1,]] - psi[pairs[2,]])^2 - (v[pairs[1,]] + v[pairs[2,]])) / 2
  max(0, mean(vals))
}

## EB mean given tau2
eb_mean <- function(psi, v, tau2) {
  w <- 1 / (v + tau2)
  sum(w * psi) / sum(w)
}

## Finite-sample conformal quantile (ceiling/floor rule from the paper).
## Returns NA when the requested rank is out of [1, J_cal].
conformal_quantile <- function(scores, p) {
  J   <- length(scores)
  idx <- if (p <= 0.5) floor(p * (J + 1)) else ceiling(p * (J + 1))
  s   <- sort(scores)
  if (idx < 1 || idx > J) return(NA_real_)
  s[idx]
}

## The minimum achievable two-sided alpha for J_cal calibration scores.
min_alpha_twosided <- function(J_cal) 2 / (J_cal + 1)

## Split conformal PI (one split).
##   Returns lo/hi = NA when alpha is too small for J_cal.
split_conformal_pi <- function(psi_T, v_T, psi_C, v_C, v_new, alpha) {
  # (i) Training
  tau2_tr   <- tau2_pairwise(psi_T, v_T)
  psi_eb_tr <- eb_mean(psi_T, v_T, tau2_tr)

  # (ii) Calibration scores
  S_cal <- (psi_C - psi_eb_tr) / sqrt(tau2_tr + v_C)

  # (iii) Quantiles
  q_lo <- conformal_quantile(S_cal, alpha / 2)
  q_hi <- conformal_quantile(S_cal, 1 - alpha / 2)

  scale <- sqrt(tau2_tr + v_new)
  lo <- if (is.na(q_lo)) NA_real_ else psi_eb_tr + q_lo * scale
  hi <- if (is.na(q_hi)) NA_real_ else psi_eb_tr + q_hi * scale

  list(lo = lo, hi = hi, width = if (any(is.na(c(lo,hi)))) NA_real_ else hi - lo,
       tau2_tr = tau2_tr, psi_eb_tr = psi_eb_tr,
       q_lo = q_lo, q_hi = q_hi, S_cal = S_cal)
}

# ── Parametric normal PI (benchmark) ─────────────────────────────────────────
parametric_pi <- function(psi, v, v_new, alpha = 0.10) {
  tau2 <- tau2_pairwise(psi, v)
  mu   <- eb_mean(psi, v, tau2)
  z    <- qnorm(1 - alpha/2)
  scale <- sqrt(tau2 + v_new)
  list(lo = mu - z * scale, hi = mu + z * scale,
       width = 2 * z * scale, tau2 = tau2, mu = mu)
}

# ── Analysis A: J=15 (tau2 > 0) ──────────────────────────────────────────────
cat("════════════════════════════════════════════\n")
cat("Analysis A: J = 15 (all functionals, tau2 > 0)\n")
cat("════════════════════════════════════════════\n\n")

psi15 <- d$psi;  v15 <- d$v
v_new <- mean(v15)   # target: new functional with average variance

# With J_train=7, J_cal=8:
#   min achievable two-sided alpha = 2/(8+1) = 2/9 ≈ 0.222  -> 77.8% PI
# With J_loo=15:
#   min achievable two-sided alpha = 2/(15+1) = 1/8 = 0.125 -> 87.5% PI
# Parametric normal PI uses 90% throughout for comparison.
alpha_norm   <- 0.10   # 90% normal PI
alpha_split  <- 2/9    # 77.8%: achievable for J_cal=8 split conformal
alpha_loo    <- 2/16   # 87.5%: achievable for J_loo=15 LOO conformal

cat(sprintf("J_cal=8 => min two-sided alpha = 2/9 = %.4f (%.1f%% PI)\n",
            alpha_split, 100*(1-alpha_split)))
cat(sprintf("J_loo=15 => min two-sided alpha = 2/16 = %.4f (%.1f%% PI)\n\n",
            alpha_loo, 100*(1-alpha_loo)))

## Parametric normal PI (90%)
pnorm15 <- parametric_pi(psi15, v15, v_new, alpha_norm)
cat(sprintf("Parametric 90%% PI (normal):        [%+7.0f, %+7.0f]  width %6.0f\n",
            pnorm15$lo, pnorm15$hi, pnorm15$width))
cat(sprintf("  tau_hat = $%.0f,  psi_EB = $%.0f\n\n",
            sqrt(pnorm15$tau2), pnorm15$mu))

## Split conformal: B random splits at the achievable 77.8% level
B <- 2000
res15 <- replicate(B, {
  idx    <- sample(15)
  J_tr   <- 7; J_cal <- 8
  T_idx  <- idx[1:J_tr];  C_idx <- idx[(J_tr+1):15]
  r      <- split_conformal_pi(psi15[T_idx], v15[T_idx],
                               psi15[C_idx], v15[C_idx],
                               v_new, alpha_split)
  c(lo=r$lo, hi=r$hi, width=r$width,
    tau2_tr=r$tau2_tr, psi_eb_tr=r$psi_eb_tr)
}, simplify = TRUE)

cat(sprintf("Split conformal %.1f%% PI (B=%d splits, J_cal=8):\n",
            100*(1-alpha_split), B))
cat(sprintf("  Median:  [%+7.0f, %+7.0f]  width %6.0f\n",
            median(res15["lo",], na.rm=TRUE),
            median(res15["hi",], na.rm=TRUE),
            median(res15["width",], na.rm=TRUE)))
cat(sprintf("  SD of endpoints: lo=%.0f  hi=%.0f\n",
            sd(res15["lo",], na.rm=TRUE), sd(res15["hi",], na.rm=TRUE)))
cat(sprintf("  Fraction of splits with finite PI: %.1f%%\n\n",
            100*mean(!is.na(res15["width",]))))

lo_med <- median(res15["lo",], na.rm=TRUE)
hi_med <- median(res15["hi",], na.rm=TRUE)
cat(sprintf("Coverage of observed psi_j by median split PI: %.0f/15 = %.1f%%\n",
            sum(psi15 >= lo_med & psi15 <= hi_med, na.rm=TRUE),
            100*mean(psi15 >= lo_med & psi15 <= hi_med, na.rm=TRUE)))
cat(sprintf("Experimental benchmark $1,794 in split PI: %s\n\n",
            ifelse(!is.na(lo_med) && 1794 >= lo_med & 1794 <= hi_med, "YES", "NO")))

# ── Analysis B: J=11 covariate-adjusted (tau2 = 0, for contrast) ─────────────
cat("════════════════════════════════════════════\n")
cat("Analysis B: J = 11 (covariate-adjusted, tau2 = 0)\n")
cat("  Note: conformal NOT valid here (tau2=0 regime uses sandwich CI)\n")
cat("════════════════════════════════════════════\n\n")

adj_idx <- !grepl("^DIM-CPS|^DIM-PSID", d$label)
psi11 <- d$psi[adj_idx];  v11 <- d$v[adj_idx]
v_new11 <- mean(v11)

pnorm11 <- parametric_pi(psi11, v11, v_new11, alpha_norm)
cat(sprintf("Parametric 90%% PI: [%+7.0f, %+7.0f]  width %6.0f\n",
            pnorm11$lo, pnorm11$hi, pnorm11$width))
cat(sprintf("  tau_hat = $%.0f,  psi_EB = $%.0f\n", sqrt(pnorm11$tau2), pnorm11$mu))
cat("  (tau2=0 => parametric PI degenerates to sandwich CI; shown for contrast)\n\n")

## J_cal=6 for J=11, so alpha_split11 = 2/7 ≈ 0.286 (71.4% PI)
alpha_split11 <- 2/7
res11 <- replicate(B, {
  idx <- sample(11)
  J_tr <- 5;  J_cal <- 6
  r <- split_conformal_pi(psi11[idx[1:J_tr]], v11[idx[1:J_tr]],
                          psi11[idx[(J_tr+1):11]], v11[idx[(J_tr+1):11]],
                          v_new11, alpha_split11)
  c(lo=r$lo, hi=r$hi, width=r$width, tau2_tr=r$tau2_tr)
}, simplify = TRUE)

cat(sprintf("  [Naive conformal %.1f%% at tau2=0 -- INVALID, for contrast]\n",
            100*(1-alpha_split11)))
cat(sprintf("  Median PI: [%+7.0f, %+7.0f]  width %6.0f  (tau_tr=%.0f)\n",
            median(res11["lo",], na.rm=TRUE), median(res11["hi",], na.rm=TRUE),
            median(res11["width",], na.rm=TRUE),
            sqrt(median(res11["tau2_tr",], na.rm=TRUE))))

# ── Analysis C: LOO diagnostic ────────────────────────────────────────────────
cat("\n════════════════════════════════════════════\n")
cat("Analysis C: Leave-one-out diagnostic (J=15)\n")
cat("════════════════════════════════════════════\n\n")

# For each j, train on J-1, compute LOO score S_j^LOO
# Then check: (a) distribution of LOO scores, (b) LOO coverage

loo <- lapply(seq_len(15), function(j) {
  T_idx <- setdiff(1:15, j)
  tau2_j  <- tau2_pairwise(psi15[T_idx], v15[T_idx])
  psi_eb_j <- eb_mean(psi15[T_idx], v15[T_idx], tau2_j)
  S_j  <- (psi15[j] - psi_eb_j) / sqrt(tau2_j + v15[j])
  data.frame(j=j, label=d$label[j], psi=psi15[j],
             tau2_tr=tau2_j, psi_eb_tr=psi_eb_j,
             S_loo=S_j,
             lo90 = psi_eb_j + qnorm(0.05) * sqrt(tau2_j + v15[j]),
             hi90 = psi_eb_j + qnorm(0.95) * sqrt(tau2_j + v15[j]))
})
loo <- do.call(rbind, loo)

cat("LOO scores S_j^LOO (should ~ N(0,1) if tau2>0 and model correct):\n")
cat(sprintf("  Mean=%.3f  SD=%.3f  (ideal: 0, 1)\n\n",
            mean(loo$S_loo), sd(loo$S_loo)))

## LOO-based PI: use LOO scores as calibration, full-data EB as center.
## Achievable level with J_loo=15: alpha_loo = 2/16 = 0.125 -> 87.5% PI
tau2_full   <- tau2_pairwise(psi15, v15)
psi_eb_full <- eb_mean(psi15, v15, tau2_full)
v_new_full  <- mean(v15)

q_lo_loo <- conformal_quantile(loo$S_loo, alpha_loo / 2)        # rank 1 of 15
q_hi_loo <- conformal_quantile(loo$S_loo, 1 - alpha_loo / 2)   # rank 15 of 15
lo_loo_pi <- psi_eb_full + q_lo_loo * sqrt(tau2_full + v_new_full)
hi_loo_pi <- psi_eb_full + q_hi_loo * sqrt(tau2_full + v_new_full)

cat(sprintf("LOO conformal %.1f%% PI for new functional:\n",
            100*(1 - alpha_loo)))
cat(sprintf("  quantiles: q_%.1f%%%% = %.3f,  q_%.1f%%%% = %.3f\n",
            100*alpha_loo/2, q_lo_loo, 100*(1-alpha_loo/2), q_hi_loo))
cat(sprintf("  PI: [%+7.0f, %+7.0f]  width %6.0f\n\n",
            lo_loo_pi, hi_loo_pi, hi_loo_pi - lo_loo_pi))

## LOO self-coverage: does normal LOO PI cover each psi_j?
loo$covered_norm <- loo$psi >= loo$lo90 & loo$psi <= loo$hi90
cat(sprintf("LOO normal PI self-coverage: %d/15 = %.0f%% (target: ~90%%)\n",
            sum(loo$covered_norm), 100*mean(loo$covered_norm)))
cat("\nFull LOO table:\n")
cat(sprintf("%-15s  %8s  %8s  %6s  %8s\n",
            "Functional", "psi_hat", "psi_EB^{-j}", "S_loo", "covered?"))
for (i in seq_len(nrow(loo))) {
  r <- loo[i,]
  cat(sprintf("%-15s  %8.0f  %8.0f  %6.2f  %s\n",
              r$label, r$psi, r$psi_eb_tr, r$S_loo,
              ifelse(r$covered_norm, "YES", "NO ")))
}

# ── Save conformal results ────────────────────────────────────────────────────
conf_out <- data.frame(
  analysis  = c("J15_normal_90",
                sprintf("J15_split_conformal_%.1f", 100*(1-alpha_split)),
                sprintf("J15_loo_conformal_%.1f",   100*(1-alpha_loo)),
                "J11_normal_90_tau0"),
  level     = c(90, 100*(1-alpha_split), 100*(1-alpha_loo), 90),
  lo        = c(pnorm15$lo,  lo_med,     lo_loo_pi,  pnorm11$lo),
  hi        = c(pnorm15$hi,  hi_med,     hi_loo_pi,  pnorm11$hi),
  width     = c(pnorm15$width, hi_med-lo_med, hi_loo_pi-lo_loo_pi, pnorm11$width),
  tau2      = c(pnorm15$tau2, pnorm15$tau2, tau2_full, pnorm11$tau2),
  psi_eb    = c(pnorm15$mu, median(res15["psi_eb_tr",],na.rm=TRUE),
                psi_eb_full, pnorm11$mu),
  valid_conformal = c(FALSE, TRUE, TRUE, FALSE)
)
write.csv(conf_out, file.path(dir_data, "lalonde_conformal.csv"), row.names = FALSE)
write.csv(loo, file.path(dir_data, "lalonde_loo.csv"), row.names = FALSE)
cat("\nSaved: lalonde_conformal.csv, lalonde_loo.csv\n")

# ── Plots ─────────────────────────────────────────────────────────────────────
library(ggplot2)

## Plot 1: calibration scores (LOO) vs N(0,1) quantiles (Q-Q plot)
loo_qq <- loo |>
  mutate(
    group = case_when(
      label %in% c("DIM-RCT","IPW-NB","OR-NB")        ~ "RCT-only",
      grepl("^DIM-", label)                             ~ "EC crude (DIM-k)",
      grepl("^OM-",  label)                             ~ "EC adjusted (OM-k)",
      grepl("^IPW-", label) & !label %in% c("IPW-NB") ~ "EC adjusted (IPW-k)"
    )
  )
palette <- c(
  "RCT-only"            = "#1f78b4",
  "EC adjusted (OM-k)"  = "#33a02c",
  "EC adjusted (IPW-k)" = "#ff7f00",
  "EC crude (DIM-k)"    = "#e31a1c"
)

p_qq <- ggplot(loo_qq |> arrange(S_loo) |> mutate(q_th = qnorm((row_number()-0.5)/n())),
               aes(x = q_th, y = S_loo, color = group, label = label)) +
  geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
  geom_hline(yintercept = c(q_lo_loo, q_hi_loo), linetype="dotted", color="grey30",
             linewidth=0.7) +
  annotate("text", x=-2.1, y=q_lo_loo+0.12,
           label=sprintf("q[%.0f%%]", 100*alpha_loo/2),
           parse=FALSE, size=2.8, hjust=0, color="grey30") +
  annotate("text", x=-2.1, y=q_hi_loo+0.12,
           label=sprintf("q[%.0f%%]", 100*(1-alpha_loo/2)),
           parse=FALSE, size=2.8, hjust=0, color="grey30") +
  geom_point(size=3) +
  ggrepel::geom_text_repel(size=2.8, max.overlaps=20, segment.size=0.2) +
  scale_color_manual(values=palette, name=NULL) +
  labs(x = "N(0,1) theoretical quantile",
       y = expression("LOO score " * S[j]^{LOO}),
       title = "Q-Q plot of LOO conformal scores (J = 15)",
       subtitle = sprintf("Dotted: %.1f%%/%.1f%% quantiles for %.1f%% LOO PI",
                          100*alpha_loo/2, 100*(1-alpha_loo/2), 100*(1-alpha_loo))) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", panel.grid.minor=element_blank(),
        plot.subtitle=element_text(size=8.5))

ggsave(file.path(dir_plots, "conformal_qq.pdf"), p_qq,
       width=6.5, height=5, device="pdf")
message("Saved conformal_qq.pdf")

## Plot 2: split-conformal PI distribution over B splits
res_df <- as.data.frame(t(res15))
p_dist <- ggplot(res_df, aes(x=lo)) +
  geom_histogram(aes(x=lo, fill="Lower endpoint"), alpha=0.6, bins=40) +
  geom_histogram(aes(x=hi, fill="Upper endpoint"), alpha=0.6, bins=40) +
  geom_vline(xintercept=1794, linetype="solid",  color="black",   linewidth=0.8) +
  geom_vline(xintercept=c(lo_med, hi_med), linetype="dashed", color="steelblue", linewidth=0.8) +
  annotate("text", x=1794+100, y=Inf, vjust=1.5, hjust=0, label="Expt.\nbenchmark",
           size=3, color="black") +
  scale_fill_manual(values=c("Lower endpoint"="#e31a1c","Upper endpoint"="#1f78b4"),
                    name=NULL) +
  labs(x="PI endpoint (USD)", y="Count",
       title=sprintf("Distribution of split conformal PI endpoints over %d random splits (J=15)", B),
       subtitle="Dashed blue: median endpoints; solid black: experimental benchmark $1,794") +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", panel.grid.minor=element_blank(),
        plot.subtitle=element_text(size=8.5))

ggsave(file.path(dir_plots, "conformal_splits.pdf"), p_dist,
       width=6.5, height=4.5, device="pdf")
message("Saved conformal_splits.pdf")

## Plot 3: PI comparison -- forest-style showing psi_j and three PI types
pi_comp <- data.frame(
  method = c(sprintf("Normal 90%% PI\n(J=15)"),
             sprintf("Split conformal %.1f%%%% PI\n(J=15, J_cal=8, median split)", 100*(1-alpha_split)),
             sprintf("LOO conformal %.1f%%%% PI\n(J=15, J_cal=15)", 100*(1-alpha_loo))),
  lo     = c(pnorm15$lo,   lo_med,     lo_loo_pi),
  hi     = c(pnorm15$hi,   hi_med,     hi_loo_pi),
  center = c(pnorm15$mu,   median(res15["psi_eb_tr",], na.rm=TRUE), psi_eb_full),
  valid  = c(FALSE, TRUE, TRUE)
)
pi_comp$method <- factor(pi_comp$method, levels=rev(pi_comp$method))

p_pi <- ggplot(pi_comp, aes(y=method)) +
  geom_rect(aes(xmin=lo, xmax=hi, ymin=as.numeric(method)-0.25, ymax=as.numeric(method)+0.25,
                fill=valid), alpha=0.25) +
  geom_segment(aes(x=lo, xend=hi, yend=method, color=valid), linewidth=1.2) +
  geom_point(aes(x=center, color=valid), size=3) +
  scale_fill_manual(values=c("TRUE"="steelblue","FALSE"="grey50"), guide="none") +
  scale_color_manual(values=c("TRUE"="steelblue","FALSE"="grey50"),
                     labels=c("TRUE"="Conformal (valid)","FALSE"="Normal (parametric)"),
                     name=NULL) +
  # Individual estimates
  geom_point(data=data.frame(psi=psi15, y=rep(0.4, 15)),
             aes(x=psi, y=y), inherit.aes=FALSE, color="grey40", size=1.5, shape=124) +
  geom_vline(xintercept=1794,  linetype="solid",  color="black", linewidth=0.7) +
  geom_vline(xintercept=684.7, linetype="dashed", color="#e31a1c", linewidth=0.7) +
  annotate("text", x=1794+80,   y=0, label="$1,794", hjust=0, size=3) +
  annotate("text", x=684.7-80,  y=0, label=expression(hat(psi)[EB]), hjust=1, size=3, color="#e31a1c") +
  scale_x_continuous(labels=scales::dollar_format(prefix="$"),
                     limits=c(-5500, 5000)) +
  labs(x="USD (1982)", y=NULL,
       title="90% prediction interval for a new functional (J=15, tau2>0)",
       subtitle="Tick marks: individual estimates. Solid: expt. benchmark. Dashed: psi_EB") +
  theme_bw(base_size=11) +
  theme(panel.grid.minor=element_blank(), plot.subtitle=element_text(size=8.5))

ggsave(file.path(dir_plots, "conformal_pi_comp.pdf"), p_pi,
       width=7, height=4, device="pdf")
message("Saved conformal_pi_comp.pdf")
