## lalonde_plots.R
## Generates figures for the LaLonde HCT case study:
##   (1) forest_plot.pdf  — 15 estimates with 95% CIs, EB lines
##   (2) eb_shrinkage.pdf — raw vs posterior mean, by functional group
##   (3) qstat.pdf        — Q statistic across three sub-analyses

library(ggplot2)
library(dplyr)

dir_data  <- "~/Dropbox/Aplicaciones/Overleaf/shrinkage/data"
dir_plots <- "~/Dropbox/Aplicaciones/Overleaf/shrinkage/plots"

d <- read.csv(file.path(dir_data, "lalonde_functionals.csv"))

# ── Metadata ─────────────────────────────────────────────────────────────────
d <- d |>
  mutate(
    group = case_when(
      label %in% c("DIM-RCT","IPW-NB","OR-NB")        ~ "RCT-only",
      grepl("^DIM-", label) & label != "DIM-RCT"       ~ "EC crude (DIM-k)",
      grepl("^OM-",  label)                             ~ "EC adjusted (OM-k)",
      grepl("^IPW-", label) & label != "IPW-NB"        ~ "EC adjusted (IPW-k)"
    ),
    group = factor(group, levels = c(
      "RCT-only", "EC adjusted (OM-k)", "EC adjusted (IPW-k)", "EC crude (DIM-k)"
    )),
    lo95  = psi - 1.96 * se,
    hi95  = psi + 1.96 * se,
    label = factor(label, levels = rev(d$label))   # bottom-to-top order
  )

palette <- c(
  "RCT-only"            = "#1f78b4",
  "EC adjusted (OM-k)"  = "#33a02c",
  "EC adjusted (IPW-k)" = "#ff7f00",
  "EC crude (DIM-k)"    = "#e31a1c"
)

benchmark  <- 1794.34
psi_eb_all <- 684.72
psi_eb_11  <- 1578.00   # from sub-analysis
psi_eb_3   <- 1796.00

# ── 1. Forest plot ────────────────────────────────────────────────────────────
p_forest <- ggplot(d, aes(x = psi, y = label, color = group)) +
  # Reference lines
  geom_vline(xintercept = benchmark,  linetype = "solid",  color = "black",  linewidth = 0.7) +
  geom_vline(xintercept = 0,          linetype = "dotted", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = psi_eb_11,  linetype = "dashed", color = "#1f78b4",linewidth = 0.8) +
  geom_vline(xintercept = psi_eb_all, linetype = "dashed", color = "#e31a1c",linewidth = 0.8) +
  # CIs and points
  geom_errorbarh(aes(xmin = lo95, xmax = hi95), height = 0.25, linewidth = 0.55) +
  geom_point(size = 2.2) +
  # Annotations
  annotate("text", x = benchmark + 60,  y = 0.55, label = "Expt. benchmark",
           hjust = 0, size = 3, color = "black") +
  annotate("text", x = psi_eb_11 + 60,  y = 1.55, label = expression(hat(psi)[EB]^{J==11}),
           hjust = 0, size = 3, color = "#1f78b4") +
  annotate("text", x = psi_eb_all - 60, y = 3.0,  label = expression(hat(psi)[EB]^{J==15}),
           hjust = 1, size = 3, color = "#e31a1c") +
  scale_color_manual(values = palette, name = NULL) +
  scale_x_continuous(labels = scales::dollar_format(prefix = "$"),
                     breaks = seq(-6000, 4000, by = 2000)) +
  labs(
    x = expression(hat(psi)[j] ~ "(earnings 1978, USD)"),
    y = NULL,
    title = "LaLonde: 15 identification functionals for the ATT",
    subtitle = expression("95% CIs via influence-function SE; vertical lines at experimental benchmark and " * hat(psi)[EB])
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position   = "bottom",
    panel.grid.minor  = element_blank(),
    plot.subtitle     = element_text(size = 8.5)
  )

ggsave(file.path(dir_plots, "forest_plot.pdf"), p_forest,
       width = 7, height = 6.5, device = "pdf")
message("Saved forest_plot.pdf")

# ── 2. EB shrinkage scatter ───────────────────────────────────────────────────
# For J=11 (adjusted only) we stored psi_EB_j per row; re-compute for J=11 sub
adj_idx <- !grepl("^DIM-CPS|^DIM-PSID", d$label)
v11 <- d$v[adj_idx]; psi11 <- d$psi[adj_idx]
w11 <- 1/v11; pbar11 <- sum(w11*psi11)/sum(w11)
Q11 <- sum((psi11-pbar11)^2/v11)
C11 <- sum(w11)-sum(w11^2)/sum(w11)
tau2_11 <- max(0, (Q11-(sum(adj_idx)-1))/C11)
web11 <- 1/(v11+tau2_11)
peb11 <- sum(web11*psi11)/sum(web11)

d$psi_eb11 <- NA_real_
d$psi_eb11[adj_idx] <- peb11 + tau2_11/(v11+tau2_11) * (psi11 - peb11)
# tau2_11 = 0, so psi_eb11 = pbar11 for all adj. functionals (uniform EB)

# Use the J=15 EB posterior means for the shrinkage scatter
d_plot2 <- d |>
  select(label, group, psi, psi_EB_j, psi_eb11) |>
  mutate(label = as.character(label))

p_shrink <- ggplot(d_plot2, aes(x = psi, y = psi_EB_j, color = group, shape = group)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey60") +
  geom_hline(yintercept = psi_eb_all, linetype = "dashed", color = "#e31a1c", linewidth = 0.6) +
  geom_hline(yintercept = benchmark,  linetype = "solid",  color = "black",  linewidth = 0.6) +
  geom_segment(aes(xend = psi, yend = psi_eb_all), linetype = "solid",
               linewidth = 0.3, alpha = 0.4) +
  geom_point(size = 3) +
  ggrepel::geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 20,
                           segment.size = 0.2) +
  scale_color_manual(values = palette, name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
  scale_x_continuous(labels = scales::dollar_format(prefix = "$"),
                     breaks = seq(-4000, 3000, by = 2000)) +
  scale_y_continuous(labels = scales::dollar_format(prefix = "$")) +
  labs(
    x = expression("Raw estimate " * hat(psi)[j]),
    y = expression("EB posterior mean " * hat(psi)[j]^{EB}),
    title = "Empirical Bayes shrinkage: raw vs posterior mean (J = 15)",
    subtitle = expression("Dashed red: " * hat(psi)[EB] == "$685" * "; solid black: expt. benchmark $1,794")
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.5))

ggsave(file.path(dir_plots, "eb_shrinkage.pdf"), p_shrink,
       width = 6.5, height = 5.5, device = "pdf")
message("Saved eb_shrinkage.pdf")

# ── 3. Q-statistic bar chart ──────────────────────────────────────────────────
qdata <- data.frame(
  analysis  = factor(
    c("RCT-only\n(J = 3)", "Covariate-adjusted\n(J = 11)", "All functionals\n(J = 15)"),
    levels = c("RCT-only\n(J = 3)", "Covariate-adjusted\n(J = 11)", "All functionals\n(J = 15)")
  ),
  Q     = c(0.00, 1.46, 99.88),
  df    = c(2, 10, 14),
  crit  = qchisq(0.95, df = c(2, 10, 14)),
  tau2  = c(0, 0, 3293908),
  psiEB = c(1796, 1578, 685)
)

p_q <- ggplot(qdata, aes(x = analysis, y = Q)) +
  geom_col(aes(fill = Q > crit), width = 0.5) +
  geom_errorbar(aes(ymin = crit, ymax = crit), width = 0.4,
                linetype = "dashed", color = "grey30", linewidth = 0.8) +
  geom_text(aes(label = sprintf("Q = %.2f\nχ²(%.0f) crit = %.1f", Q, df, crit)),
            vjust = -0.4, size = 3.0) +
  geom_text(aes(y = -15,
                label = sprintf("hat(psi)[EB] == '$%.0f'", psiEB)),
            parse = TRUE, size = 3.2, color = "grey20") +
  scale_fill_manual(values = c("FALSE" = "#1f78b4", "TRUE" = "#e31a1c"),
                    labels = c("Fail to reject H0: tau^2 = 0",
                               "Reject H0: tau^2 > 0"),
                    name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.12, 0.15))) +
  labs(
    x = NULL,
    y = "Cochran Q statistic",
    title = "Heterogeneity test across three sub-analyses",
    subtitle = "Dashed bars: chi-squared 95% critical value (df = J - 1)"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.5))

ggsave(file.path(dir_plots, "qstat.pdf"), p_q,
       width = 6.5, height = 4.5, device = "pdf")
message("Saved qstat.pdf")

# ── 4. Variance decomposition (v_j vs tau2) ───────────────────────────────────
tau2_all <- unique(d$tau2_hat)
d_var <- d |>
  mutate(
    sigma2 = v,
    signal_frac = tau2_all / (v + tau2_all),
    label = as.character(label)
  )

p_var <- ggplot(d_var, aes(x = reorder(label, se), y = se, fill = group)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = sqrt(tau2_all), linetype = "dashed",
             color = "#e31a1c", linewidth = 0.8) +
  annotate("text", x = 1, y = sqrt(tau2_all) + 30,
           label = expression(hat(tau) == "$1815"), hjust = 0, size = 3, color = "#e31a1c") +
  scale_fill_manual(values = palette, name = NULL) +
  scale_y_continuous(labels = scales::dollar_format(prefix = "$")) +
  labs(
    x = NULL,
    y = expression("Standard error " * sqrt(v[j])),
    title = "Standard errors of individual functionals vs heterogeneity scale",
    subtitle = expression("Dashed: " * hat(tau) * " (scale of between-functional spread, J = 15)")
  ) +
  coord_flip() +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(size = 8.5))

ggsave(file.path(dir_plots, "vardecomp.pdf"), p_var,
       width = 6.5, height = 5.5, device = "pdf")
message("Saved vardecomp.pdf")
