# =============================================================================
# 07_FAAH_NAAH_sex_comparison.R
# Confronto maschi vs femmine Tg2576 per gli enzimi di degradazione della PEA
#
# Razionale:
#   I risultati su NOR, spine e 3-NT mostrano una risposta attenuata alla PEA
#   nelle femmine Tg2576 rispetto ai maschi. La PEA endogena ed esogena viene
#   degradata da due enzimi principali: FAAH e NAAA. Se la loro espressione
#   è maggiore nelle femmine, la PEA somministrata verrebbe metabolizzata
#   più rapidamente, riducendone l'efficacia terapeutica.
#
# Ipotesi:
#   FAAH e NAAA sono più espressi nelle femmine Tg2576 rispetto ai maschi.
#
# Confronti (ciascuno a 2 gruppi):
#   1. FAAH mRNA    — Tg_male_veh vs Tg_female_veh
#   2. FAAH protein — Tg_male_veh vs Tg_female_veh
#   3. NAAH mRNA    — Tg_male_veh vs Tg_female_veh
#
# Scelta test: Shapiro-Wilk → parametrico o non parametrico
# Nessuna correzione per confronti multipli: ciascun enzima/readout è un
# esperimento indipendente con ipotesi propria (α = 0.05 per confronto)
# Effect size: Cohen's d o Cliff's delta
# =============================================================================

# --- Pacchetti ---------------------------------------------------------------
required_pkgs <- c("readr", "tidyr", "dplyr", "effsize")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(readr)
library(tidyr)
library(dplyr)
library(effsize)

# --- Funzione di analisi per un enzima ---------------------------------------
analyze_enzyme <- function(filename, enzyme_label, alpha = 0.05) {

  df <- read_csv(file.path("DATASET", filename), show_col_types = FALSE)

  long <- df |>
    pivot_longer(cols = everything(), names_to = "group", values_to = "value") |>
    filter(!is.na(value))

  males   <- long$value[long$group == "Tg_male_veh"]
  females <- long$value[long$group == "Tg_female_veh"]

  cat("=======================================================================\n")
  cat(sprintf("  %s — Tg MASCHI vs Tg FEMMINE\n", toupper(enzyme_label)))
  cat(sprintf("  Ipotesi: espressione maggiore nelle femmine\n"))
  cat(sprintf("  Test indipendente (α = %.2f)\n", alpha))
  cat("=======================================================================\n\n")

  # --- Statistiche descrittive -----------------------------------------------
  cat("--- Statistiche descrittive ---\n\n")
  desc <- long |>
    group_by(group) |>
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

  # Fold change (femmine / maschi)
  fc <- mean(females) / mean(males)
  cat(sprintf("\n  Fold change (femmine / maschi) = %.2f\n", fc))

  # --- Test di normalità -----------------------------------------------------
  cat("\n--- Test di normalità (Shapiro-Wilk) ---\n\n")
  norm_m <- shapiro.test(males)
  norm_f <- shapiro.test(females)
  cat(sprintf("  Maschi:  W = %.3f, p = %.4f  %s\n",
              norm_m$statistic, norm_m$p.value,
              ifelse(norm_m$p.value > 0.05, "(normale)", "(NON normale)")))
  cat(sprintf("  Femmine: W = %.3f, p = %.4f  %s\n",
              norm_f$statistic, norm_f$p.value,
              ifelse(norm_f$p.value > 0.05, "(normale)", "(NON normale)")))

  all_normal <- (norm_m$p.value > 0.05) & (norm_f$p.value > 0.05)

  # --- Test statistico -------------------------------------------------------
  if (all_normal) {
    cat("\n→ Test parametrico: Welch t-test + Cohen's d\n\n")

    # Direzione: femmine vs maschi (t positivo → femmine > maschi)
    tt <- t.test(females, males)
    cat(sprintf("Welch t-test: t = %.3f, df = %.1f, p = %.4f\n",
                tt$statistic, tt$parameter, tt$p.value))
    cat(sprintf("  Mean maschi  = %.6f (SD = %.6f, n = %d)\n",
                mean(males), sd(males), length(males)))
    cat(sprintf("  Mean femmine = %.6f (SD = %.6f, n = %d)\n",
                mean(females), sd(females), length(females)))
    cat(sprintf("  Differenza media (femmine - maschi) = %.6f\n",
                mean(females) - mean(males)))
    cat(sprintf("  IC 95%%: [%.6f, %.6f]\n",
                tt$conf.int[1], tt$conf.int[2]))

    cat(sprintf("\n  p < %.2f? %s\n",
                alpha,
                ifelse(tt$p.value < alpha, "SI", "NO")))

    d <- cohen.d(females, males)
    cat(sprintf("\n  Cohen's d = %.3f  [IC 95%%: %.3f, %.3f]\n",
                d$estimate, d$conf.int[1], d$conf.int[2]))
    cat(sprintf("  Interpretazione: %s\n", as.character(d$magnitude)))
    cat("  (d positivo → femmine > maschi)\n")

    return(invisible(list(
      test = "t-test", stat = tt$statistic, p = tt$p.value,
      es = d$estimate, es_ci = d$conf.int,
      es_mag = as.character(d$magnitude), es_type = "Cohen's d",
      fc = fc, desc = desc
    )))

  } else {
    cat("\n→ Test non parametrico: Mann-Whitney U + Cliff's delta\n\n")

    wt <- wilcox.test(females, males, conf.int = TRUE)
    cat(sprintf("Mann-Whitney U: W = %.1f, p = %.4f\n",
                wt$statistic, wt$p.value))
    cat(sprintf("  Median maschi  = %.6f (IQR = %.6f, n = %d)\n",
                median(males), IQR(males), length(males)))
    cat(sprintf("  Median femmine = %.6f (IQR = %.6f, n = %d)\n",
                median(females), IQR(females), length(females)))
    cat(sprintf("  Hodges-Lehmann estimate = %.6f\n", wt$estimate))
    cat(sprintf("  IC 95%%: [%.6f, %.6f]\n",
                wt$conf.int[1], wt$conf.int[2]))

    cat(sprintf("\n  p < %.2f? %s\n",
                alpha,
                ifelse(wt$p.value < alpha, "SI", "NO")))

    cd <- cliff.delta(females, males)
    cat(sprintf("\n  Cliff's delta = %.3f  [IC 95%%: %.3f, %.3f]\n",
                cd$estimate, cd$conf.int[1], cd$conf.int[2]))
    cat(sprintf("  Interpretazione: %s\n", as.character(cd$magnitude)))
    cat("  (delta positivo → femmine > maschi)\n")

    return(invisible(list(
      test = "Mann-Whitney", stat = wt$statistic, p = wt$p.value,
      es = cd$estimate, es_ci = cd$conf.int,
      es_mag = as.character(cd$magnitude), es_type = "Cliff's d",
      fc = fc, desc = desc
    )))
  }
}

# =============================================================================
# ESECUZIONE
# =============================================================================

if (!dir.exists("OUTPUT")) dir.create("OUTPUT")

output_file <- "OUTPUT/07_FAAH_NAAH_sex_comparison.txt"
sink(output_file, split = TRUE)
cat("Analisi eseguita il:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# --- FAAH mRNA ---------------------------------------------------------------
r_faah_mrna <- analyze_enzyme("faah_genico.csv", "FAAH mRNA")

# --- FAAH proteina -----------------------------------------------------------
r_faah_prot <- analyze_enzyme("faah_proteica.csv", "FAAH proteina")

# --- NAAH mRNA ---------------------------------------------------------------
r_naah_mrna <- analyze_enzyme("naah_genico.csv", "NAAH mRNA")

# --- Tabella riassuntiva -----------------------------------------------------
cat("\n=======================================================================\n")
cat("  TABELLA RIASSUNTIVA — ENZIMI DI DEGRADAZIONE DELLA PEA\n")
cat("  Confronto: Tg femmine vs Tg maschi (vehicle)\n")
cat("=======================================================================\n\n")

results <- list(
  list(name = "FAAH mRNA",     res = r_faah_mrna),
  list(name = "FAAH proteina", res = r_faah_prot),
  list(name = "NAAH mRNA",     res = r_naah_mrna)
)

cat(sprintf("%-16s %-12s %8s %10s %10s %6s %s\n",
            "Enzima", "Test", "Stat", "p", "ES", "FC", "p < 0.05"))
cat(strrep("-", 82), "\n")

for (r in results) {
  cat(sprintf("%-16s %-12s %8.3f %10.4f %10.3f %6.2f %s\n",
              r$name,
              r$res$es_type,
              r$res$stat,
              r$res$p,
              r$res$es,
              r$res$fc,
              ifelse(r$res$p < 0.05, "SI *", "NO")))
}

cat("\nFC = fold change (femmine / maschi). FC > 1 = espressione maggiore nelle femmine.\n")
cat("ES positivo = femmine > maschi.\n")

sink()

# --- Salva oggetti -----------------------------------------------------------
save(r_faah_mrna, r_faah_prot, r_naah_mrna,
     file = "OUTPUT/07_FAAH_NAAH_sex_comparison.RData")

cat("\nFile salvati in OUTPUT/:\n")
cat(" -", basename(output_file), "\n")
cat(" - 07_FAAH_NAAH_sex_comparison.RData\n")
