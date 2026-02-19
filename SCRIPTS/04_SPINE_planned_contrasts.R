# =============================================================================
# 04_SPINE_planned_contrasts.R
# Confronti pianificati su Spine Density (Apicale e Basale) — Solo Tg2576
#
# Ipotesi a priori (basata su Tortolani et al., 2025):
#   La PEA ripristina la densità delle spine dendritiche nei Tg2576.
#   Ci aspettiamo un effetto positivo sia nei maschi che nelle femmine.
#
# Confronti per ciascun compartimento (apicale, basale):
#   1. Tg_male_veh   vs Tg_male_PEA   (replicazione dato pubblicato)
#   2. Tg_female_veh vs Tg_female_PEA  (estensione alle femmine)
#
# Nota: diversi gruppi non superano Shapiro-Wilk → test non parametrici
#   Mann-Whitney U + Cliff's delta come effect size
# Correzione: Bonferroni per 2 confronti per compartimento (α = 0.025)
# =============================================================================

# --- Pacchetti ---------------------------------------------------------------
required_pkgs <- c("readr", "tidyr", "dplyr", "effsize")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(readr)
library(tidyr)
library(dplyr)
library(effsize)

# --- Funzione di analisi per un compartimento --------------------------------
analyze_spine <- function(filename, compartment_label) {

  df <- read_csv(file.path("DATASET", filename), show_col_types = FALSE)

  long <- df |>
    pivot_longer(cols = everything(), names_to = "group", values_to = "value") |>
    filter(!is.na(value)) |>
    filter(grepl("^Tg", group)) |>
    mutate(
      sex       = ifelse(grepl("female", group), "female", "male"),
      treatment = ifelse(grepl("PEA", group), "PEA", "veh")
    )

  cat("=======================================================================\n")
  cat(sprintf("  SPINE DENDRITICHE %s — CONFRONTI PIANIFICATI (Solo Tg2576)\n",
              toupper(compartment_label)))
  cat("  Ipotesi a priori: PEA ripristina la densità in entrambi i sessi\n")
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
    summarise(
      n    = n(),
      SW_W = shapiro.test(value)$statistic,
      SW_p = shapiro.test(value)$p.value,
      .groups = "drop"
    )
  print(norm)
  all_normal <- all(norm$SW_p > 0.05)
  cat(sprintf("\nTutti i gruppi normali (p > 0.05): %s\n", all_normal))

  if (!all_normal) {
    cat("→ Almeno un gruppo non normale: si usano test non parametrici\n")
    cat("  (Mann-Whitney U + Cliff's delta)\n")
  }

  # --- Kruskal-Wallis omnibus (contesto) -------------------------------------
  cat("\n=======================================================================\n")
  cat("  KRUSKAL-WALLIS (analisi omnibus di contesto)\n")
  cat("=======================================================================\n\n")

  kw <- kruskal.test(value ~ interaction(sex, treatment), data = long)
  cat(sprintf("Kruskal-Wallis: chi-sq = %.3f, df = %d, p = %.6f\n",
              kw$statistic, kw$parameter, kw$p.value))

  # --- CONFRONTO 1: Maschi ---------------------------------------------------
  cat("\n=======================================================================\n")
  cat(sprintf("  CONFRONTO 1: Tg MASCHI %s — veh vs PEA\n", compartment_label))
  cat("  (Replicazione del dato pubblicato, Tortolani et al. 2025)\n")
  cat("=======================================================================\n\n")

  m_veh <- long$value[long$sex == "male" & long$treatment == "veh"]
  m_pea <- long$value[long$sex == "male" & long$treatment == "PEA"]

  # Mann-Whitney U
  w_male <- wilcox.test(m_pea, m_veh, conf.int = TRUE)
  cat(sprintf("Mann-Whitney U: W = %.1f, p = %.4f\n",
              w_male$statistic, w_male$p.value))
  cat(sprintf("  Median veh = %.3f (IQR = %.3f, n = %d)\n",
              median(m_veh), IQR(m_veh), length(m_veh)))
  cat(sprintf("  Median PEA = %.3f (IQR = %.3f, n = %d)\n",
              median(m_pea), IQR(m_pea), length(m_pea)))
  cat(sprintf("  Hodges-Lehmann estimate = %.3f\n", w_male$estimate))
  cat(sprintf("  IC 95%%: [%.3f, %.3f]\n",
              w_male$conf.int[1], w_male$conf.int[2]))

  cat(sprintf("\n  p < 0.025 (Bonferroni)? %s\n",
              ifelse(w_male$p.value < 0.025, "SI", "NO")))

  # Cliff's delta
  cd_male <- cliff.delta(m_pea, m_veh)
  cat(sprintf("\n  Cliff's delta = %.3f  [IC 95%%: %.3f, %.3f]\n",
              cd_male$estimate, cd_male$conf.int[1], cd_male$conf.int[2]))
  cat(sprintf("  Interpretazione: %s\n", as.character(cd_male$magnitude)))

  # --- CONFRONTO 2: Femmine -------------------------------------------------
  cat("\n=======================================================================\n")
  cat(sprintf("  CONFRONTO 2: Tg FEMMINE %s — veh vs PEA\n", compartment_label))
  cat("  (Estensione alle femmine — dato nuovo)\n")
  cat("=======================================================================\n\n")

  f_veh <- long$value[long$sex == "female" & long$treatment == "veh"]
  f_pea <- long$value[long$sex == "female" & long$treatment == "PEA"]

  # Mann-Whitney U
  w_female <- wilcox.test(f_pea, f_veh, conf.int = TRUE)
  cat(sprintf("Mann-Whitney U: W = %.1f, p = %.4f\n",
              w_female$statistic, w_female$p.value))
  cat(sprintf("  Median veh = %.3f (IQR = %.3f, n = %d)\n",
              median(f_veh), IQR(f_veh), length(f_veh)))
  cat(sprintf("  Median PEA = %.3f (IQR = %.3f, n = %d)\n",
              median(f_pea), IQR(f_pea), length(f_pea)))
  cat(sprintf("  Hodges-Lehmann estimate = %.3f\n", w_female$estimate))
  cat(sprintf("  IC 95%%: [%.3f, %.3f]\n",
              w_female$conf.int[1], w_female$conf.int[2]))

  cat(sprintf("\n  p < 0.025 (Bonferroni)? %s\n",
              ifelse(w_female$p.value < 0.025, "SI", "NO")))

  # Cliff's delta
  cd_female <- cliff.delta(f_pea, f_veh)
  cat(sprintf("\n  Cliff's delta = %.3f  [IC 95%%: %.3f, %.3f]\n",
              cd_female$estimate, cd_female$conf.int[1], cd_female$conf.int[2]))
  cat(sprintf("  Interpretazione: %s\n", as.character(cd_female$magnitude)))

  # --- Riepilogo -------------------------------------------------------------
  cat("\n=======================================================================\n")
  cat(sprintf("  RIEPILOGO — %s\n", toupper(compartment_label)))
  cat("=======================================================================\n\n")

  cat(sprintf("%-20s %8s %8s %12s %10s %s\n",
              "Confronto", "W", "p", "Cliff's d", "Magnitude", "Bonf.(p<.025)"))
  cat(strrep("-", 78), "\n")
  cat(sprintf("%-20s %8.1f %8.4f %12.3f %10s %s\n",
              "Maschi (PEA vs veh)",
              w_male$statistic, w_male$p.value,
              cd_male$estimate, as.character(cd_male$magnitude),
              ifelse(w_male$p.value < 0.025, "SI *", "NO")))
  cat(sprintf("%-20s %8.1f %8.4f %12.3f %10s %s\n",
              "Femmine (PEA vs veh)",
              w_female$statistic, w_female$p.value,
              cd_female$estimate, as.character(cd_female$magnitude),
              ifelse(w_female$p.value < 0.025, "SI *", "NO")))

  cat("\n-----------------------------------------------------------------------\n")
  cat("Interpretazione Cliff's delta:\n")
  cat("  |d| < 0.147 = trascurabile, 0.147-0.330 = piccolo,\n")
  cat("  0.330-0.474 = medio, > 0.474 = grande\n")
  cat("-----------------------------------------------------------------------\n\n")

  return(invisible(list(
    desc = desc, norm = norm,
    w_male = w_male, w_female = w_female,
    cd_male = cd_male, cd_female = cd_female
  )))
}

# =============================================================================
# ESECUZIONE
# =============================================================================

if (!dir.exists("OUTPUT")) dir.create("OUTPUT")

output_file <- "OUTPUT/04_SPINE_planned_contrasts.txt"
sink(output_file, split = TRUE)
cat("Analisi eseguita il:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- Spine Apicali -----------------------------------------------------------
r_apical <- analyze_spine("apical_spine.csv", "apicali (CA1)")

# --- Spine Basali ------------------------------------------------------------
r_basal  <- analyze_spine("basal_spine.csv",  "basali (CA1)")

sink()

# --- Salva oggetti -----------------------------------------------------------
save(r_apical, r_basal, file = "OUTPUT/04_SPINE_planned_contrasts.RData")

cat("\nFile salvati in OUTPUT/:\n")
cat(" -", basename(output_file), "\n")
cat(" - 04_SPINE_planned_contrasts.RData\n")
