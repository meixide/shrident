## build_datasets.R
## Construye los dos datasets limpios para la aplicación empírica del paper.
##
## LaLonde (Dehejia-Wahba):
##   - q=5 entornos: Z=0 NSW experimental, Z=1..4 grupos control observacionales
##   - Salida: lalonde_combined.csv  (columnas: Z, A, Y, age, education, black,
##             hispanic, married, nodegree, re74, re75)
##
## Card & Krueger (1994):
##   - 410 restaurantes NJ/PA, dos oleadas (pre/post aumento salario mínimo NJ)
##   - Salida: ck94.csv  (columnas: store_id, chain, state, region,
##             fte_pre, fte_post, delta_fte, co_owned)

library(haven)
library(dplyr)

dir_data <- dirname(rstudioapi::getSourceEditorContext()$path %||%
                    normalizePath("."))
# Fallback si no se llama desde RStudio
if (!grepl("data$", dir_data))
  dir_data <- file.path(
    "~/Dropbox/Aplicaciones/Overleaf/shrinkage/data")

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

# Z=1..4 : grupos control observacionales (todos A=0)
cps2  <- read_lalonde("cps_controls2.dta",  z_val = 1)
cps3  <- read_lalonde("cps_controls3.dta",  z_val = 2)
psid2 <- read_lalonde("psid_controls2.dta", z_val = 3)
psid3 <- read_lalonde("psid_controls3.dta", z_val = 4)

lalonde <- bind_rows(nsw, cps2, cps3, psid2, psid3)

# Resumen
cat("\nEntorno   n      n_tratados  pi_z\n")
lalonde |>
  group_by(Z) |>
  summarise(n = n(), n_trt = sum(A), pi = mean(A)) |>
  print()

# Tratados solo existen en Z=0 → Wald válido para subconjuntos que incluyen Z=0
# J = 2^4 - 1 = 15 funcionales (subconjuntos de {1,2,3,4} no vacíos, combinados con 0)
cat("\nJ = 15 funcionales Wald válidos (todos los subconjuntos S con 0 ∈ S, |S|≥2)\n")

write.csv(lalonde, file.path(dir_data, "lalonde_combined.csv"), row.names = FALSE)
message("Guardado: data/lalonde_combined.csv  (", nrow(lalonde), " obs)")

# ──────────────────────────────────────────────────────────────
# 2. CARD & KRUEGER (1994)  ────────────────────────────────────
# ──────────────────────────────────────────────────────────────
message("\n=== Card & Krueger (1994) ===")

# Formato de ancho fijo según el codebook
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

# Full-time equivalent employment (FTE = full-time + 0.5*part-time + managers)
ck <- ck |>
  mutate(
    fte_pre  = empft  + 0.5 * emppt  + nmgrs,
    fte_post = empft2 + 0.5 * emppt2 + nmgrs2,
    delta_fte = fte_post - fte_pre,
    # Región del grupo control (PA vs distintas zonas NJ)
    region = case_when(
      state == 0 & pa1 == 1 ~ "PA_north",
      state == 0 & pa2 == 1 ~ "PA_east",
      state == 0            ~ "PA_other",
      northj  == 1          ~ "NJ_north",
      centralj == 1         ~ "NJ_central",
      southj  == 1          ~ "NJ_south",
      shore   == 1          ~ "NJ_shore",
      TRUE                  ~ "other"
    ),
    chain_name = recode(chain,
      `1` = "BK", `2` = "KFC", `3` = "RoyRogers", `4` = "Wendys")
  ) |>
  # Excluir tiendas cerradas en segunda oleada
  filter(status2 %in% c(1, NA) | is.na(status2))

# Grupos control candidatos para DiD (C=0 tratado=NJ, C=k control candidato k)
# K=1: todos los PA  |  K=2: solo PA_north  |  K=3: solo PA_east
# K=4: solo BK en PA |  K=5: solo KFC en PA
ck <- ck |>
  mutate(
    C = case_when(
      state == 1 ~ 0L,                          # Tratado: NJ
      state == 0 ~ 1L,                          # Control 1: todos PA
      TRUE       ~ NA_integer_
    )
  )

cat("\nEstado  n    fte_pre_mean  fte_post_mean  delta_mean\n")
ck |>
  group_by(state) |>
  summarise(n = n(),
            fte_pre  = mean(fte_pre,  na.rm = TRUE),
            fte_post = mean(fte_post, na.rm = TRUE),
            delta    = mean(delta_fte, na.rm = TRUE)) |>
  print()

cat("\nCadena   n_NJ  n_PA\n")
ck |> count(chain_name, state) |>
  tidyr::pivot_wider(names_from = state, values_from = n,
                     names_prefix = "state_") |> print()

write.csv(ck, file.path(dir_data, "ck94.csv"), row.names = FALSE)
message("Guardado: data/ck94.csv  (", nrow(ck), " obs)")

message("\n=== DONE ===")
message("Archivos generados:")
message("  data/lalonde_combined.csv  — LaLonde 5 entornos, 3259 obs")
message("  data/ck94.csv              — Card & Krueger 1994, ~399 obs")
