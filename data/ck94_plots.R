## ck94_plots.R
## Generates figures for the Card & Krueger (1994) DiD case study:
##   (1) ck94_forest.pdf      — 25 estimates with 95% CIs, EB line
##   (2) ck94_qstat.pdf       — Q statistic across DIM, OM, IPW sub-analyses
##   (3) ck94_shrinkage.pdf   — raw vs posterior mean, by method
##   (4) ck94_vardecomp.pdf   — SE bar chart vs tau_hat
##   (5) ck94_conformal_qq.pdf — LOO conformal score Q-Q plot

library(ggplot2)
library(dplyr)

dir_data  <- "~/Dropbox/Aplicaciones/Overleaf/shrinkage/data"
dir_plots <- "~/Dropbox/Aplicaciones/Overleaf/shrinkage/plots"

d <- read.csv(file.path(dir_data, "ck94_functionals.csv")) |>
  mutate(
    method = factor(method, levels = c("DIM", "OM", "IPW")),
    label  = factor(label, levels = rev(label))
  )

palette <- c(
  "DIM" = "#e31a1c",
  "OM"  = "#33a02c",
  "IPW" = "#ff7f00"
)

benchmark <- 2.75       # DIM-PA-all (original CK94)
psi_eb    <- unique(d$psiEB)
tau_hat   <- sqrt(unique(d$tau2_hat))

# ── 1. Forest plot ─────────────────────────────────────────────────────────────
p_forest <- ggplot(d, aes(x = psi, y = label, color = method)) +
  geom_vline(xintercept = 0,         linetype = "dotted", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = benchmark, linetype = "solid",  color = "black",  linewidth = 0.7) +
  geom_vline(xintercept = psi_eb,    linetype = "dashed", color = "#1f78b4",linewidth = 0.8) +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95), height = 0.25, linewidth = 0.55) +
  geom_point(size = 2.2) +
  annotate("text", x = benchmark + 0.15, y = 0.6,
           label = "DIM-PA-all\n(CK94 benchmark)",
           hjust = 0, size = 2.8, color = "black") +
  annotate("text", x = psi_eb - 0.15, y = 3.5,
           label = expression(hat(psi)[EB]^{J==25}),
           hjust = 1, size = 3, color = "#1f78b4") +
  scale_color_manual(values = palette, name = "Method") +
  scale_x_continuous(
    breaks = seq(-4, 10, by = 2),
    labels = function(x) paste0(x, " FTE")
  ) +
  labs(
    x = expression(hat(psi)[j] ~ "(change in FTE employment)"),
    y = NULL,
    title = "Card & Krueger (1994): 25 DiD functionals for the ATT",
    subtitle = expression("95% CIs via HC variance; vertical lines at CK94 benchmark and " * hat(psi)[EB])
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.subtitle    = element_text(size = 8.5)
  )

ggsave(file.path(dir_plots, "ck94_forest.pdf"), p_forest,
       width = 7, height = 7.5, device = "pdf")
message("Saved ck94_forest.pdf")

# ── 2. Q-statistic bar chart ──────────────────────────────────────────────────
cochran_Q <- function(psi, v) {
  w <- 1/v; pbar <- sum(w*psi)/sum(w)
  sum(w*(psi-pbar)^2)
}
mom_tau2 <- function(psi, v) {
  w <- 1/v; pbar <- sum(w*psi)/sum(w)
  Q <- sum(w*(psi-pbar)^2)
  C <- sum(w) - sum(w^2)/sum(w)
  max(0, (Q-(length(psi)-1))/C)
}
eb_mean <- function(psi, v, tau2) {
  web <- 1/(v+tau2); sum(web*psi)/sum(web)
}

psi_all <- d$psi; v_all <- d$v
psi_dim <- d$psi[d$method=="DIM"]; v_dim <- d$v[d$method=="DIM"]
psi_om  <- d$psi[d$method=="OM"]; v_om <- d$v[d$method=="OM"]
psi_ipw <- d$psi[d$method=="IPW"]; v_ipw <- d$v[d$method=="IPW"]

qdata <- data.frame(
  analysis = factor(
    c("DIM only\n(J = 9)", "OM only\n(J = 8)", "IPW only\n(J = 8)", "All methods\n(J = 25)"),
    levels = c("DIM only\n(J = 9)", "OM only\n(J = 8)", "IPW only\n(J = 8)", "All methods\n(J = 25)")
  ),
  Q    = c(cochran_Q(psi_dim,v_dim), cochran_Q(psi_om,v_om),
           cochran_Q(psi_ipw,v_ipw), cochran_Q(psi_all,v_all)),
  df   = c(8, 7, 7, 24),
  tau2 = c(mom_tau2(psi_dim,v_dim), mom_tau2(psi_om,v_om),
           mom_tau2(psi_ipw,v_ipw), mom_tau2(psi_all,v_all)),
  psiEB = c(eb_mean(psi_dim,v_dim,mom_tau2(psi_dim,v_dim)),
            eb_mean(psi_om,v_om,mom_tau2(psi_om,v_om)),
            eb_mean(psi_ipw,v_ipw,mom_tau2(psi_ipw,v_ipw)),
            eb_mean(psi_all,v_all,mom_tau2(psi_all,v_all)))
)
qdata$crit <- qchisq(0.95, df = qdata$df)

p_q <- ggplot(qdata, aes(x = analysis, y = Q)) +
  geom_col(aes(fill = Q > crit), width = 0.55) +
  geom_errorbar(aes(ymin = crit, ymax = crit), width = 0.45,
                linetype = "dashed", color = "grey30", linewidth = 0.8) +
  geom_text(aes(label = sprintf("Q = %.1f\nχ²(%d) = %.1f", Q, df, crit)),
            vjust = -0.3, size = 2.8) +
  geom_text(aes(y = -25,
                label = sprintf("hat(psi)[EB] == '%.2f FTE'", psiEB)),
            parse = TRUE, size = 3.0, color = "grey20") +
  scale_fill_manual(
    values = c("FALSE" = "#1f78b4", "TRUE" = "#e31a1c"),
    labels = c("Fail to reject H₀: τ² = 0",
               "Reject H₀: τ² > 0"),
    name = NULL
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
  labs(
    x = NULL,
    y = "Cochran Q statistic",
    title = "Heterogeneity test across three sub-analyses",
    subtitle = "Dashed bars: chi-squared 95% critical value (df = J - 1)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.5))

ggsave(file.path(dir_plots, "ck94_qstat.pdf"), p_q,
       width = 7, height = 4.5, device = "pdf")
message("Saved ck94_qstat.pdf")

# ── 3. EB shrinkage scatter ───────────────────────────────────────────────────
p_shrink <- ggplot(d, aes(x = psi, y = psi_EB_j, color = method, shape = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
  geom_hline(yintercept = psi_eb,    linetype = "dashed", color = "#1f78b4", linewidth = 0.6) +
  geom_hline(yintercept = benchmark, linetype = "solid",  color = "black",  linewidth = 0.6) +
  geom_segment(aes(xend = psi, yend = psi_eb), linetype = "solid",
               linewidth = 0.3, alpha = 0.35) +
  geom_point(size = 2.8) +
  ggrepel::geom_text_repel(aes(label = label), size = 2.4, max.overlaps = 25,
                            segment.size = 0.2) +
  scale_color_manual(values = palette, name = "Method") +
  scale_shape_manual(values = c(16, 17, 15), name = "Method") +
  scale_x_continuous(labels = function(x) paste0(x, " FTE"),
                     breaks = seq(-4, 10, by = 2)) +
  scale_y_continuous(labels = function(x) paste0(x, " FTE")) +
  labs(
    x = expression("Raw estimate " * hat(psi)[j]),
    y = expression("EB posterior mean " * hat(psi)[j]^{EB}),
    title = "Empirical Bayes shrinkage: raw vs posterior mean (J = 25)",
    subtitle = sprintf("Dashed blue: psi_EB = %.2f FTE; solid black: CK94 benchmark 2.75 FTE", psi_eb)
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.5))

ggsave(file.path(dir_plots, "ck94_shrinkage.pdf"), p_shrink,
       width = 6.5, height = 5.5, device = "pdf")
message("Saved ck94_shrinkage.pdf")

# ── 4. Variance decomposition ─────────────────────────────────────────────────
p_var <- ggplot(d |> mutate(label = as.character(label)),
               aes(x = reorder(label, se), y = se, fill = method)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = tau_hat, linetype = "dashed",
             color = "#1f78b4", linewidth = 0.8) +
  annotate("text", x = 1.5, y = tau_hat + 0.08,
           label = sprintf("hat(tau) == '%.2f FTE'", tau_hat),
           parse = TRUE, hjust = 0, size = 3, color = "#1f78b4") +
  scale_fill_manual(values = palette, name = "Method") +
  scale_y_continuous(labels = function(x) paste0(x, " FTE")) +
  labs(
    x = NULL,
    y = expression("Standard error " * sqrt(v[j])),
    title = "Standard errors vs heterogeneity scale (J = 25)",
    subtitle = expression("Dashed: " * hat(tau) * " (between-control-group spread)")
  ) +
  coord_flip() +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.5))

ggsave(file.path(dir_plots, "ck94_vardecomp.pdf"), p_var,
       width = 6.5, height = 6, device = "pdf")
message("Saved ck94_vardecomp.pdf")

# ── 5. LOO conformal Q-Q ─────────────────────────────────────────────────────
J <- nrow(d)
# Theoretical quantiles from uniform order statistics
probs <- (seq_len(J) - 0.5) / J
s_sorted <- sort(d$S_loo)

# Standard normal reference
qq_df <- data.frame(
  theoretical = qnorm(probs),
  observed    = s_sorted,
  method      = d$method[order(d$S_loo)]
)

p_qq <- ggplot(qq_df, aes(x = theoretical, y = observed, color = method)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5) +
  scale_color_manual(values = palette, name = "Method") +
  labs(
    x = "Standard normal quantile",
    y = expression("LOO conformal score " * S[j]^{LOO}),
    title = "Q-Q plot of LOO conformal scores (J = 25)",
    subtitle = "Dashed: reference N(0,1); heavy tails indicate tau2 > 0"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.5))

ggsave(file.path(dir_plots, "ck94_conformal_qq.pdf"), p_qq,
       width = 5.5, height = 4.5, device = "pdf")
message("Saved ck94_conformal_qq.pdf")
