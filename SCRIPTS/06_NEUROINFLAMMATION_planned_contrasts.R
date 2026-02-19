# =============================================================================
# 06_NEUROINFLAMMATION_planned_contrasts.R
# Confronti pianificati su marcatori neuroinfiammatori — Solo Tg2576
#
# Marcatori:
#   GFAP  — astrogliosi
#   IBA1  — microgliosi
#   iNOS  — enzima produttore di NO (legato a stress nitrosativo / 3-NT)
#
# Ipotesi a priori (basata su Tortolani et al., 2025):
#   La PEA riduce la neuroinfiammazione nei Tg2576.
#   Ci aspettiamo una riduzione in entrambi i sessi.
#
# Confronti per ciascun marcatore:
#   1. Tg_male_veh   vs Tg_male_PEA   (replicazione dato pubblicato)
#   2. Tg_female_veh vs Tg_female_PEA  (estensione alle femmine)
#
# Scelta test: Shapiro-Wilk → parametrico (Welch t + Cohen's d)
#                              o non parametrico (Mann-Whitney U + Cliff's delta)
# Correzione: Bonferroni per 2 confronti per marcatore (α = 0.025)
# =============================================================================

# --- Pacchetti ---------------------------------------------------------------
required_pkgs <- c("readr", "tidyr", "dplyr", "effsize")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(readr)
library(tidyr)
library(dplyr)
library(effsize)

# --- Funzione di analisi generica -------------------------------------------
analyze_marker <- function(filename, marker_label, direction = "lower_is_better") {
  # direction: "lower_is_better" → veh > PEA = PEA riduce il marcatore (atteso)
  #            "higher_is_better" → PEA > veh

  df <- read_csv(file.path("DATASET", filename), show_col_types = FALSE)

  long <- df |>
    pivot_longer(cols = everything(), names_to = "group", values_to = "value") |>
    filter(!is.na(value)) |>
    mutate(
      sex       = ifelse(grepl("female", group), "female", "male"),
      treatment = ifelse(grepl("PEA", group), "PEA", "veh")
    )

  cat("=======================================================================\n")
  cat(sprintf("  %s — CONFRONTI PIANIFICATI (Solo Tg2576)\n", toupper(marker_label)))
  cat("  Ipotesi a priori: PEA riduce il marcatore in entrambi i sessi\n")
  cat("  Correzione Bonferroni per 2 confronti (α = 0.025)\n")
  cat("=======================================================================\n\n")

  # --- Statistiche descrittive -----------------------------------------------
  cat("--- Statistiche descrittive ---\n\n")
  desc <- long |>
    group_by(sex, treatment) |>
    summarise(
      n      = n(),
      mean   = mean(value),
      sd     = sd(value),
      median = median(value),
      min    = min(value),
      max    = max(value),
      .groups = "drop"
    )
  print(desc)

  # --- Test di normalità -----------------------------------------------------
  cat("\n--- Test di normalità (Shapiro-Wilk) ---\n\n")
  norm <- long |>
    group_by(sex, treatment) |>
    filter(n() >= 3) |>
    summarise(
      n    = n(),
      SW_W = shapiro.test(value)$statistic,
      SW_p = shapiro.test(value)$p.value,
      .groups = "drop"
    )
  print(norm)
  all_normal <- all(norm$SW_p > 0.05)
  cat(sprintf("\nTutti i gruppi normali (p > 0.05): %s\n", all_normal))

  use_parametric <- all_normal
  if (use_parametric) {
    cat("→ Test parametrici: Welch t-test + Cohen's d\n")
  } else {
    cat("→ Test non parametrici: Mann-Whitney U + Cliff's delta\n")
  }

  # --- Nota su n piccoli -----------------------------------------------------
  n_min <- min(desc$n)
  if (n_min <= 5) {
    cat(sprintf("\n⚠ ATTENZIONE: n minimo = %d. Potenza statistica limitata.\n", n_min))
  }

  # --- Funzione interna per un confronto ------------------------------------
  run_contrast <- function(sex_label, sex_val, contrast_title) {
    cat("\n=======================================================================\n")
    cat(sprintf("  CONFRONTO: Tg %s — veh vs PEA\n", toupper(sex_label)))
    cat(sprintf("  %s\n", contrast_title))
    cat("=======================================================================\n\n")

    grp_veh <- long$value[long$sex == sex_val & long$treatment == "veh"]
    grp_pea <- long$value[long$sex == sex_val & long$treatment == "PEA"]

    n_veh <- length(grp_veh)
    n_pea <- length(grp_pea)

    if (use_parametric) {
      # Welch t-test (veh vs PEA: t positivo → veh > PEA → PEA riduce)
      tt <- t.test(grp_veh, grp_pea)
      cat(sprintf("Welch t-test: t = %.3f, df = %.1f, p = %.4f\n",
                  tt$statistic, tt$parameter, tt$p.value))
      cat(sprintf("  Mean veh = %.6f (SD = %.6f, n = %d)\n",
                  mean(grp_veh), sd(grp_veh), n_veh))
      cat(sprintf("  Mean PEA = %.6f (SD = %.6f, n = %d)\n",
                  mean(grp_pea), sd(grp_pea), n_pea))
      cat(sprintf("  Differenza media (veh - PEA) = %.6f\n",
                  mean(grp_veh) - mean(grp_pea)))
      cat(sprintf("  IC 95%%: [%.6f, %.6f]\n",
                  tt$conf.int[1], tt$conf.int[2]))

      cat(sprintf("\n  p < 0.025 (Bonferroni)? %s\n",
                  ifelse(tt$p.value < 0.025, "SI", "NO")))

      d <- cohen.d(grp_veh, grp_pea)
      cat(sprintf("\n  Cohen's d = %.3f  [IC 95%%: %.3f, %.3f]\n",
                  d$estimate, d$conf.int[1], d$conf.int[2]))
      cat(sprintf("  Interpretazione: %s\n", as.character(d$magnitude)))

      return(list(test = "t-test", stat = tt$statistic, p = tt$p.value,
                  es = d$estimate, es_ci = d$conf.int,
                  es_mag = as.character(d$magnitude), es_type = "Cohen's d"))

    } else {
      # Mann-Whitney U
      wt <- wilcox.test(grp_veh, grp_pea, conf.int = TRUE)
      cat(sprintf("Mann-Whitney U: W = %.1f, p = %.4f\n",
                  wt$statistic, wt$p.value))
      cat(sprintf("  Median veh = %.6f (IQR = %.6f, n = %d)\n",
                  median(grp_veh), IQR(grp_veh), n_veh))
      cat(sprintf("  Median PEA = %.6f (IQR = %.6f, n = %d)\n",
                  median(grp_pea), IQR(grp_pea), n_pea))
      cat(sprintf("  Hodges-Lehmann estimate = %.6f\n", wt$estimate))
      cat(sprintf("  IC 95%%: [%.6f, %.6f]\n",
                  wt$conf.int[1], wt$conf.int[2]))

      cat(sprintf("\n  p < 0.025 (Bonferroni)? %s\n",
                  ifelse(wt$p.value < 0.025, "SI", "NO")))

      cd <- cliff.delta(grp_veh, grp_pea)
      cat(sprintf("\n  Cliff's delta = %.3f  [IC 95%%: %.3f, %.3f]\n",
                  cd$estimate, cd$conf.int[1], cd$conf.int[2]))
      cat(sprintf("  Interpretazione: %s\n", as.character(cd$magnitude)))

      return(list(test = "Mann-Whitney", stat = wt$statistic, p = wt$p.value,
                  es = cd$estimate, es_ci = cd$conf.int,
                  es_mag = as.character(cd$magnitude), es_type = "Cliff's d"))
    }
  }

  # --- Esegui confronti -----------------------------------------------------
  res_male   <- run_contrast("maschi",  "male",
                             "(Replicazione del dato pubblicato, Tortolani et al. 2025)")
  res_female <- run_contrast("femmine", "female",
                             "(Estensione alle femmine — dato nuovo)")

  # --- Riepilogo -------------------------------------------------------------
  cat("\n=======================================================================\n")
  cat(sprintf("  RIEPILOGO — %s\n", toupper(marker_label)))
  cat("=======================================================================\n\n")

  es_label <- res_male$es_type
  cat(sprintf("%-20s %8s %8s %12s %10s %s\n",
              "Confronto", "Stat", "p", es_label, "Magnitude", "Bonf.(p<.025)"))
  cat(strrep("-", 78), "\n")
  cat(sprintf("%-20s %8.3f %8.4f %12.3f %10s %s\n",
              "Maschi (veh vs PEA)",
              res_male$stat, res_male$p,
              res_male$es, res_male$es_mag,
              ifelse(res_male$p < 0.025, "SI *", "NO")))
  cat(sprintf("%-20s %8.3f %8.4f %12.3f %10s %s\n",
              "Femmine (veh vs PEA)",
              res_female$stat, res_female$p,
              res_female$es, res_female$es_mag,
              ifelse(res_female$p < 0.025, "SI *", "NO")))
  cat("\n")

  return(invisible(list(
    desc = desc, norm = norm,
    male = res_male, female = res_female,
    parametric = use_parametric
  )))
}

# =============================================================================
# ESECUZIONE
# =============================================================================

if (!dir.exists("OUTPUT")) dir.create("OUTPUT")

output_file <- "OUTPUT/06_NEUROINFLAMMATION_planned_contrasts.txt"
sink(output_file, split = TRUE)
cat("Analisi eseguita il:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- GFAP (astrogliosi) -----------------------------------------------------
r_gfap <- analyze_marker("gfap.csv", "GFAP (astrogliosi)")

# --- IBA1 (microgliosi) -----------------------------------------------------
r_iba1 <- analyze_marker("iba1.csv", "IBA1 (microgliosi)")

# --- iNOS (produzione NO → stress nitrosativo) ------------------------------
r_inos <- analyze_marker("inos.csv", "iNOS (produzione NO)")

# --- Tabella riassuntiva globale ---------------------------------------------
cat("=======================================================================\n")
cat("  TABELLA RIASSUNTIVA — MARCATORI NEUROINFIAMMATORI\n")
cat("=======================================================================\n\n")

markers <- list(
  list(name = "GFAP",  res = r_gfap),
  list(name = "IBA1",  res = r_iba1),
  list(name = "iNOS",  res = r_inos)
)

cat(sprintf("%-8s %-12s | %-10s %8s %10s | %-10s %8s %10s\n",
            "", "",
            "MASCHI", "p", "ES",
            "FEMMINE", "p", "ES"))
cat(strrep("-", 85), "\n")

for (m in markers) {
  cat(sprintf("%-8s %-12s | %-10s %8.4f %10.3f | %-10s %8.4f %10.3f\n",
              m$name,
              m$res$male$es_type,
              m$res$male$es_mag,
              m$res$male$p,
              m$res$male$es,
              m$res$female$es_mag,
              m$res$female$p,
              m$res$female$es))
}
cat("\n")

sink()

# --- Salva oggetti -----------------------------------------------------------
save(r_gfap, r_iba1, r_inos,
     file = "OUTPUT/06_NEUROINFLAMMATION_planned_contrasts.RData")

cat("\nFile salvati in OUTPUT/:\n")
cat(" -", basename(output_file), "\n")
cat(" - 06_NEUROINFLAMMATION_planned_contrasts.RData\n")
