## ck_forest.R
## Forest plot of the 25 Card-Krueger estimates with jackknife intervals.
##
##   Rscript apps/ck_forest.R results/ck94_intervals_jackknife_B2000_s42
##
## Reads estimates.csv and pooled.csv of a run of apps/ck_intervals.py --jackknife and
## writes ck94_forest.pdf into the same directory. Vertical lines: solid, the estimate
## against all of Pennsylvania (DIM-PA-all); dashed, the precision-weighted pooled
## estimate; dotted, the convex combination.
suppressMessages({library(ggplot2); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
run <- args[1]
d <- read.csv(file.path(run, "estimates.csv")) |>
  mutate(method = factor(method, levels = c("DIM", "OM", "IPW")),
         label = factor(label, levels = rev(label)), lo = psi - 1.96 * se, hi = psi + 1.96 * se)
p <- read.csv(file.path(run, "pooled.csv")) |> filter(set == "All")
bench <- d$psi[d$label == "DIM-PA-all"]
g <- ggplot(d, aes(x = psi, y = label, shape = method)) +
  geom_vline(xintercept = bench, linetype = "solid", linewidth = 0.5) +
  geom_vline(xintercept = p$psi_eb_zero, linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = p$psi_convex_zero, linetype = "dotted", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.45, colour = "grey35") +
  geom_point(size = 2) +
  scale_shape_manual(values = c(DIM = 16, OM = 17, IPW = 15), name = NULL) +
  scale_x_continuous(breaks = seq(-8, 16, by = 4)) +
  labs(x = "Change in full-time-equivalent employment", y = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(run, "ck94_forest.pdf"), g, width = 7, height = 7.5, device = "pdf")
cat("written", file.path(run, "ck94_forest.pdf"), "\n")
