## build_datasets.R
## Builds the two clean data sets used in the applications of the article.
## Run from the repository root: Rscript data/build_datasets.R
##
## LaLonde (Dehejia-Wahba):
##   - 5 samples: Z=0 NSW experiment, Z=1..4 observational comparison samples
##   - Output: lalonde_combined.csv  (columns: Z, A, Y, age, education, black,
##             hispanic, married, nodegree, re74, re75)
##
## Card & Krueger (1994):
##   - 410 NJ/PA restaurants, two waves (before/after the NJ minimum-wage increase)
##   - Output: ck94.csv  (columns of public.dat as named in the codebook, plus
##             fte_pre, fte_post, delta_fte)

library(haven)
library(dplyr)

dir_data <- "data"  # relative to the repository root

# ──────────────────────────────────────────────────────────────
# 1. LALONDE  ──────────────────────────────────────────────────
# ──────────────────────────────────────────────────────────────
message("=== LaLonde / Dehejia-Wahba ===")

read_lalonde <- function(f, z_val) {
  d <- read_dta(file.path(dir_data, "lalonde", f))
  d |>
    mutate(Z = z_val) |>
    select(Z, A = treat, Y = re78,
           age, education, black, hispanic, married, nodegree, re74, re75)
}

# Z=0 : NSW experimental (treated + control from DW subsample)
nsw <- read_lalonde("nsw_dw.dta", z_val = 0)

# Z=1..4 : observational comparison samples (all A=0)
cps2  <- read_lalonde("cps_controls2.dta",  z_val = 1)
cps3  <- read_lalonde("cps_controls3.dta",  z_val = 2)
psid2 <- read_lalonde("psid_controls2.dta", z_val = 3)
psid3 <- read_lalonde("psid_controls3.dta", z_val = 4)

lalonde <- bind_rows(nsw, cps2, cps3, psid2, psid3)

# Summary
cat("\nSample    n      n_treated   pi_z\n")
lalonde |>
  group_by(Z) |>
  summarise(n = n(), n_trt = sum(A), pi = mean(A)) |>
  print()

# Treated units exist only in Z=0, so each comparison must include Z=0
# J = 2^4 - 1 = 15 functionals (non-empty subsets of {1,2,3,4}, combined with 0)
cat("\nJ = 15 Wald functionals (all subsets S with 0 in S, |S| >= 2)\n")

write.csv(lalonde, file.path(dir_data, "lalonde_combined.csv"), row.names = FALSE)
message("Saved: data/lalonde_combined.csv  (", nrow(lalonde), " obs)")

# ──────────────────────────────────────────────────────────────
# 2. CARD & KRUEGER (1994)  ────────────────────────────────────
# ──────────────────────────────────────────────────────────────
message("\n=== Card & Krueger (1994) ===")

# Fixed-width format, as in the codebook
col_spec <- list(
  sheet    = c(1,  3),
  chain    = c(5,  5),
  co_owned = c(7,  7),
  state    = c(9,  9),
  southj   = c(11, 11),
  centralj = c(13, 13),
  northj   = c(15, 15),
  pa1      = c(17, 17),
  pa2      = c(19, 19),
  shore    = c(21, 21),
  empft    = c(26, 30),
  emppt    = c(32, 36),
  nmgrs    = c(38, 42),
  wage_st  = c(44, 48),
  status2  = c(110, 110),
  empft2   = c(122, 126),
  emppt2   = c(128, 132),
  nmgrs2   = c(134, 138),
  wage_st2 = c(140, 144)
)

raw <- readLines(file.path(dir_data, "card_krueger", "public.dat"))
raw <- raw[nchar(trimws(raw)) > 0]

parse_row <- function(line) {
  vals <- lapply(col_spec, function(pos) {
    s <- substr(line, pos[1], pos[2])
    suppressWarnings(as.numeric(trimws(s)))
  })
  as.data.frame(vals)
}

ck <- do.call(rbind, lapply(raw, parse_row))

# Full-time equivalent employment (FTE = full-time + 0.5*part-time + managers).
# All 410 restaurants are kept; the analysis drops those without both FTE values.
ck <- ck |>
  select(sheet, chain, co_owned, state, southj, centralj, northj, pa1, pa2, shore,
         empft, emppt, nmgrs, status2, empft2, emppt2, nmgrs2) |>
  mutate(
    fte_pre   = empft  + 0.5 * emppt  + nmgrs,
    fte_post  = empft2 + 0.5 * emppt2 + nmgrs2,
    delta_fte = fte_post - fte_pre
  )

cat("\nState   n    fte_pre_mean  fte_post_mean  delta_mean\n")
ck |>
  group_by(state) |>
  summarise(n = n(),
            fte_pre  = mean(fte_pre,  na.rm = TRUE),
            fte_post = mean(fte_post, na.rm = TRUE),
            delta    = mean(delta_fte, na.rm = TRUE)) |>
  print()

write.csv(ck, file.path(dir_data, "ck94.csv"), row.names = FALSE)
message("Saved: data/ck94.csv  (", nrow(ck), " obs)")

message("\n=== DONE ===")
message("Files written:")
message("  data/lalonde_combined.csv  — LaLonde, 5 samples, 3624 rows")
message("  data/ck94.csv              — Card & Krueger 1994, 410 restaurants")
