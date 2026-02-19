# =============================================================================
# 02_statistics.R
# Analisi statistiche per tutti i dataset
# Approccio ibrido: test normalità (Shapiro-Wilk) → parametrico o non parametrico
#
# Strategia:
#   2 gruppi : t-test (se normale) / Mann-Whitney U (se non normale)
#   4 gruppi : ANOVA 2 vie + Tukey (se normale) / Kruskal-Wallis + Dunn (se non normale)
#   8 gruppi : ANOVA 2 vie per genotipo + Tukey / Kruskal-Wallis + Dunn
# =============================================================================

# --- Pacchetti ---------------------------------------------------------------
required_pkgs <- c("readr", "tidyr", "dplyr", "rstatix", "dunn.test")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(readr)
library(tidyr)
library(dplyr)
library(rstatix)
library(dunn.test)

# --- Funzioni di supporto ----------------------------------------------------

# Converte un dataset wide in long con colonna 'group'
to_long <- function(df) {
  df |>
    pivot_longer(cols = everything(), names_to = "group", values_to = "value") |>
    filter(!is.na(value))
}

# Test di normalità Shapiro-Wilk per ogni gruppo (n >= 3)
test_normality <- function(long_df) {
  long_df |>
    group_by(group) |>
    filter(n() >= 3) |>
    summarise(
      n       = n(),
      SW_W    = shapiro.test(value)$statistic,
      SW_p    = shapiro.test(value)$p.value,
      normal  = shapiro.test(value)$p.value > 0.05,
      .groups = "drop"
    )
}

# Decide se usare test parametrico (tutti i gruppi normali) o non parametrico
is_parametric <- function(normality_df) {
  all(normality_df$normal)
}

# Stampa risultati formattati
print_section <- function(title) {
  cat("\n", strrep("=", 70), "\n")
  cat(" ", title, "\n")
  cat(strrep("=", 70), "\n")
}

# --- Analisi 2 gruppi --------------------------------------------------------
# Confronto Tg_male_veh vs Tg_female_veh
# Dataset: faah_genico, faah_proteica, naah_genico

analyze_2groups <- function(filename, label) {
  print_section(label)
  df   <- read_csv(file.path("DATASET", filename), show_col_types = FALSE)
  long <- to_long(df)
  norm <- test_normality(long)

  cat("\nTest di normalità (Shapiro-Wilk):\n")
  print(norm)

  if (is_parametric(norm)) {
    cat("\n→ Distribuzione normale in tutti i gruppi: t-test indipendente\n")
    x <- long$value[long$group == "Tg_male_veh"]
    y <- long$value[long$group == "Tg_female_veh"]
    res <- t.test(x, y)
    cat(sprintf("  t = %.3f, df = %.1f, p = %.4f\n", res$statistic, res$parameter, res$p.value))
    cat(sprintf("  Mean Tg_male_veh = %.4f | Mean Tg_female_veh = %.4f\n",
                mean(x), mean(y)))
  } else {
    cat("\n→ Distribuzione non normale in almeno un gruppo: Mann-Whitney U\n")
    x <- long$value[long$group == "Tg_male_veh"]
    y <- long$value[long$group == "Tg_female_veh"]
    res <- wilcox.test(x, y)
    cat(sprintf("  W = %.1f, p = %.4f\n", res$statistic, res$p.value))
    cat(sprintf("  Median Tg_male_veh = %.4f | Median Tg_female_veh = %.4f\n",
                median(x), median(y)))
  }

  # Statistiche descrittive
  cat("\nStatistiche descrittive:\n")
  long |>
    group_by(group) |>
    summarise(n = n(), mean = mean(value), sd = sd(value),
              median = median(value), .groups = "drop") |>
    print()

  return(invisible(list(normality = norm, data = long)))
}

# --- Analisi 4 gruppi --------------------------------------------------------
# Fattori: Sex (male/female) × Treatment (veh/PEA)
# Dataset: 3nt, iba1, gfap, inos

analyze_4groups <- function(filename, label) {
  print_section(label)
  df   <- read_csv(file.path("DATASET", filename), show_col_types = FALSE)
  long <- to_long(df)
  norm <- test_normality(long)

  cat("\nTest di normalità (Shapiro-Wilk):\n")
  print(norm)

  # Estrai fattori da nome gruppo
  long <- long |>
    mutate(
      sex       = ifelse(grepl("male", group), "male", "female"),
      sex       = ifelse(grepl("female", group), "female", sex),
      treatment = ifelse(grepl("PEA", group), "PEA", "veh")
    )

  if (is_parametric(norm)) {
    cat("\n→ Distribuzione normale: ANOVA 2 vie (Sex × Treatment) + Tukey\n")
    model <- aov(value ~ sex * treatment, data = long)
    cat("\nANOVA:\n")
    print(summary(model))
    cat("\nPost-hoc Tukey HSD:\n")
    print(TukeyHSD(model))
  } else {
    cat("\n→ Distribuzione non normale: Kruskal-Wallis + post-hoc Dunn\n")
    kw <- kruskal.test(value ~ group, data = long)
    cat(sprintf("  Kruskal-Wallis: chi-sq = %.3f, df = %d, p = %.4f\n",
                kw$statistic, kw$parameter, kw$p.value))
    cat("\nPost-hoc Dunn (Bonferroni):\n")
    dunn.test(long$value, long$group, method = "bonferroni", kw = FALSE)
  }

  # Statistiche descrittive
  cat("\nStatistiche descrittive:\n")
  long |>
    group_by(group) |>
    summarise(n = n(), mean = mean(value), sd = sd(value),
              median = median(value), .groups = "drop") |>
    print()

  return(invisible(list(normality = norm, data = long)))
}

# --- Analisi 8 gruppi --------------------------------------------------------
# ANOVA 2 vie separata per genotipo (WT e Tg)
# Fattori: Sex × Treatment
# Dataset: nor_test, apical_spine, basal_spine

analyze_8groups <- function(filename, label) {
  print_section(label)
  df   <- read_csv(file.path("DATASET", filename), show_col_types = FALSE)
  long <- to_long(df)
  norm <- test_normality(long)

  cat("\nTest di normalità (Shapiro-Wilk):\n")
  print(norm)

  # Estrai fattori
  long <- long |>
    mutate(
      genotype  = ifelse(grepl("WT", group), "WT", "Tg"),
      sex       = ifelse(grepl("female", group), "female", "male"),
      treatment = ifelse(grepl("PEA", group), "PEA", "veh")
    )

  for (geno in c("WT", "Tg")) {
    cat(sprintf("\n--- Genotipo: %s ---\n", geno))
    sub <- long |> filter(genotype == geno)
    sub_norm <- test_normality(sub |> select(group, value))

    if (is_parametric(sub_norm)) {
      cat("→ Distribuzione normale: ANOVA 2 vie (Sex × Treatment) + Tukey\n")
      model <- aov(value ~ sex * treatment, data = sub)
      print(summary(model))
      cat("\nPost-hoc Tukey HSD:\n")
      print(TukeyHSD(model))
    } else {
      cat("→ Distribuzione non normale: Kruskal-Wallis + post-hoc Dunn\n")
      kw <- kruskal.test(value ~ group, data = sub)
      cat(sprintf("  Kruskal-Wallis: chi-sq = %.3f, df = %d, p = %.4f\n",
                  kw$statistic, kw$parameter, kw$p.value))
      cat("\nPost-hoc Dunn (Bonferroni):\n")
      dunn.test(sub$value, sub$group, method = "bonferroni", kw = FALSE)
    }
  }

  # Statistiche descrittive
  cat("\nStatistiche descrittive:\n")
  long |>
    group_by(group) |>
    summarise(n = n(), mean = mean(value), sd = sd(value),
              median = median(value), .groups = "drop") |>
    print()

  return(invisible(list(normality = norm, data = long)))
}

# =============================================================================
# ESECUZIONE ANALISI
# =============================================================================

# --- Crea cartella OUTPUT se non esiste --------------------------------------
if (!dir.exists("OUTPUT")) dir.create("OUTPUT")

# --- Avvia cattura output testuale -------------------------------------------
output_txt <- file.path("OUTPUT", "risultati_statistiche.txt")
sink(output_txt, split = TRUE)   # split=TRUE: stampa anche in console
cat("Analisi eseguita il:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# --- 2 gruppi: enzimi di degradazione della PEA ------------------------------
r_faah_genico   <- analyze_2groups("faah_genico.csv",   "FAAH mRNA (Tg_male vs Tg_female)")
r_faah_proteica <- analyze_2groups("faah_proteica.csv", "FAAH proteina (Tg_male vs Tg_female)")
r_naah_genico   <- analyze_2groups("naah_genico.csv",   "NAAH mRNA (Tg_male vs Tg_female)")

# --- 4 gruppi: marcatori infiammatori/ossidativi -----------------------------
r_3nt  <- analyze_4groups("3nt.csv",   "3-Nitrotyrosine (Sex × Treatment)")
r_iba1 <- analyze_4groups("iba1.csv",  "IBA1 (Sex × Treatment)")
r_gfap <- analyze_4groups("gfap.csv",  "GFAP (Sex × Treatment)")
r_inos <- analyze_4groups("inos.csv",  "iNOS (Sex × Treatment)")

# --- 8 gruppi: endpoint comportamentali e strutturali ------------------------
r_nor          <- analyze_8groups("nor_test.csv",     "NOR Test (Genotipo × Sex × Treatment)")
r_apical_spine <- analyze_8groups("apical_spine.csv", "Apical Spine Density (Genotipo × Sex × Treatment)")
r_basal_spine  <- analyze_8groups("basal_spine.csv",  "Basal Spine Density (Genotipo × Sex × Treatment)")

cat("\n", strrep("=", 70), "\n")
cat("  ANALISI COMPLETATA\n")
cat(strrep("=", 70), "\n")

# --- Chiudi cattura output testuale ------------------------------------------
sink()

# --- Salva oggetti R in file .RData ------------------------------------------
output_rdata <- file.path("OUTPUT", "risultati_statistiche.RData")
save(r_faah_genico, r_faah_proteica, r_naah_genico,
     r_3nt, r_iba1, r_gfap, r_inos,
     r_nor, r_apical_spine, r_basal_spine,
     file = output_rdata)

cat("\nFile salvati in OUTPUT/:\n")
cat(" -", basename(output_txt),   "\n")
cat(" -", basename(output_rdata), "\n")
