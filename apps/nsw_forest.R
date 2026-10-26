## nsw_forest.R
## Forest plot of the article's NSW estimators with person-level sandwich intervals.
##
##   Rscript apps/nsw_forest.R results/nsw_if_hapc_K5_s20261012
##
## Reads assembled/article/estimators.csv and pooled.csv (apps/nsw_article_numbers.py)
## and writes assembled/article/nsw_forest.pdf. Vertical lines: solid, experimental
## difference in means; dashed, precision-weighted pooled estimate of the
## covariate-adjusted set; dotted, its convex combination.
suppressMessages({library(ggplot2); library(dplyr)})
run <- commandArgs(trailingOnly = TRUE)[1]
dir <- file.path(run, "assembled", "article")
d <- read.csv(file.path(dir, "estimators.csv")) |>
  mutate(kind = factor(ifelse(grepl("^DIM", label), "Difference in means",
                       ifelse(grepl("^DRext", label), "Doubly robust, external regression", "Doubly robust, pooled or trial regression")),
                       levels = c("Difference in means", "Doubly robust, pooled or trial regression", "Doubly robust, external regression")),
         label = factor(label, levels = rev(label)), lo = psi - 1.96 * se, hi = psi + 1.96 * se)
p <- read.csv(file.path(dir, "pooled.csv")) |> filter(set == "onestep_xfit", analysis == "Covariate-adjusted")
g <- ggplot(d, aes(x = psi, y = label, shape = kind)) +
  geom_vline(xintercept = d$psi[d$label == "DIM-RCT"], linetype = "solid", linewidth = 0.5) +
  geom_vline(xintercept = p$psi_eb_zero, linetype = "dashed", linewidth = 0.5) +
  geom_vline(xintercept = p$psi_convex_zero, linetype = "dotted", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.25, linewidth = 0.45, colour = "grey35") +
  geom_point(size = 2) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  labs(x = "Effect on 1978 earnings of participants (dollars)", y = NULL) +
  guides(shape = guide_legend(nrow = 3)) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
ggsave(file.path(dir, "nsw_forest.pdf"), g, width = 7, height = 6, device = "pdf")
cat("written", file.path(dir, "nsw_forest.pdf"), "\n")
