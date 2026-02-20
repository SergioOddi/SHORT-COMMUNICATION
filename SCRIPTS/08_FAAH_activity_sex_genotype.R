# =============================================================================
# 08_FAAH_activity_sex_genotype.R
# Attività enzimatica FAAH — confronto Sesso × Genotipo (solo veicolo)
#
# Razionale:
#   L'espressione di FAAH (mRNA e proteina) è significativamente maggiore
#   nelle femmine Tg2576 rispetto ai maschi (sezione 3.4). Per verificare
#   se la sovraespressione si traduce in una maggiore capacità catabolica,
#   misuriamo l'attività enzimatica FAAH nell'ippocampo.
#   Includiamo anche i WT per valutare se il dimorfismo sessuale è intrinseco
#   (presente anche nei WT) o specifico del genotipo transgenico.
#
# Design:
#   2 fattori between-subject: Genotipo (Tg, WT) × Sesso (maschio, femmina)
#   Variabile dipendente: attività enzimatica FAAH (arbitrary units)
#   Solo animali vehicle-treated
#
# Confronti pianificati:
#   1. Tg maschi vs Tg femmine   (effetto sesso nel transgenico)
#   2. WT maschi vs WT femmine   (effetto sesso nel wild-type)
#   Correzione Bonferroni per 2 confronti (α = 0.025)
#
# Scelta test: Shapiro-Wilk → parametrico o non parametrico
# Effect size: Cohen's d o Cliff's delta
# =============================================================================

# --- Pacchetti ---------------------------------------------------------------
required_pkgs <- c("readr", "tidyr", "dplyr", "effsize")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs, repos = "https://cloud.r-project.org")

library(readr)
library(tidyr)
library(dplyr)
library(effsize)

# --- Caricamento dati --------------------------------------------------------
df <- read_csv("DATASET/activity.csv", show_col_types = FALSE)

long <- df |>
  pivot_longer(cols = everything(), names_to = "group", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(
    group     = trimws(group),
    genotype  = ifelse(grepl("^Tg", group), "Tg", "WT"),
    sex       = ifelse(grepl("female", group), "female", "male")
  )

# --- Output ------------------------------------------------------------------
if (!dir.exists("OUTPUT")) dir.create("OUTPUT")

output_file <- "OUTPUT/08_FAAH_activity_sex_genotype.txt"
sink(output_file, split = TRUE)

cat("=======================================================================\n")
cat("  ATTIVITÀ ENZIMATICA FAAH — Sesso × Genotipo (solo veicolo)\n")
cat("  Ipotesi: attività FAAH maggiore nelle femmine (coerente con\n")
cat("  sovraespressione mRNA/proteina osservata nella sezione 3.4)\n")
cat("  Correzione Bonferroni per 2 confronti pianificati (α = 0.025)\n")
cat("=======================================================================\n")
cat("Analisi eseguita il:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# --- Statistiche descrittive -------------------------------------------------
cat("\n--- Statistiche descrittive ---\n\n")
desc <- long |>
  group_by(genotype, sex) |>
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

# Fold change femmine / maschi per genotipo
fc_tg <- mean(long$value[long$genotype == "Tg" & long$sex == "female"]) /
         mean(long$value[long$genotype == "Tg" & long$sex == "male"])
fc_wt <- mean(long$value[long$genotype == "WT" & long$sex == "female"]) /
         mean(long$value[long$genotype == "WT" & long$sex == "male"])
cat(sprintf("\n  Fold change (femmine / maschi) — Tg = %.2f\n", fc_tg))
cat(sprintf("  Fold change (femmine / maschi) — WT = %.2f\n", fc_wt))

# --- Test di normalità -------------------------------------------------------
cat("\n--- Test di normalità (Shapiro-Wilk) ---\n\n")
norm <- long |>
  group_by(genotype, sex) |>
  summarise(
    n    = n(),
    SW_W = shapiro.test(value)$statistic,
    SW_p = shapiro.test(value)$p.value,
    .groups = "drop"
  )
print(norm)
all_normal <- all(norm$SW_p > 0.05)
cat(sprintf("\nTutti i gruppi normali (p > 0.05): %s\n", all_normal))

# --- ANOVA 2 vie (Sesso × Genotipo) ------------------------------------------
cat("\n=======================================================================\n")
cat("  ANOVA 2 VIE (Sex × Genotype) — analisi omnibus\n")
cat("=======================================================================\n\n")

model <- aov(value ~ sex * genotype, data = long)
print(summary(model))

ss <- summary(model)[[1]]
ss_total <- sum(ss[["Sum Sq"]])
cat("\nEta-squared (η²):\n")
cat(sprintf("  Sex:              η² = %.4f\n", ss[["Sum Sq"]][1] / ss_total))
cat(sprintf("  Genotype:         η² = %.4f\n", ss[["Sum Sq"]][2] / ss_total))
cat(sprintf("  Sex × Genotype:   η² = %.4f\n", ss[["Sum Sq"]][3] / ss_total))

# --- CONFRONTO 1: Tg maschi vs Tg femmine ------------------------------------
cat("\n=======================================================================\n")
cat("  CONFRONTO 1: Tg2576 — MASCHI vs FEMMINE\n")
cat("  (Conferma funzionale della sovraespressione FAAH)\n")
cat("=======================================================================\n\n")

tg_m <- long$value[long$genotype == "Tg" & long$sex == "male"]
tg_f <- long$value[long$genotype == "Tg" & long$sex == "female"]

# Check normalità per scelta test
norm_tg_m <- shapiro.test(tg_m)
norm_tg_f <- shapiro.test(tg_f)
tg_normal <- (norm_tg_m$p.value > 0.05) & (norm_tg_f$p.value > 0.05)

if (tg_normal) {
  cat("→ Test parametrico: Welch t-test + Cohen's d\n\n")
  t_tg <- t.test(tg_f, tg_m)
  cat(sprintf("Welch t-test: t = %.3f, df = %.1f, p = %.4f\n",
              t_tg$statistic, t_tg$parameter, t_tg$p.value))
  cat(sprintf("  Mean maschi  = %.4f (SD = %.4f, n = %d)\n",
              mean(tg_m), sd(tg_m), length(tg_m)))
  cat(sprintf("  Mean femmine = %.4f (SD = %.4f, n = %d)\n",
              mean(tg_f), sd(tg_f), length(tg_f)))
  cat(sprintf("  Differenza media (femmine - maschi) = %.4f\n",
              mean(tg_f) - mean(tg_m)))
  cat(sprintf("  IC 95%%: [%.4f, %.4f]\n",
              t_tg$conf.int[1], t_tg$conf.int[2]))
  cat(sprintf("\n  p < 0.025 (Bonferroni)? %s\n",
              ifelse(t_tg$p.value < 0.025, "SI", "NO")))

  d_tg <- cohen.d(tg_f, tg_m)
  cat(sprintf("\n  Cohen's d = %.3f  [IC 95%%: %.3f, %.3f]\n",
              d_tg$estimate, d_tg$conf.int[1], d_tg$conf.int[2]))
  cat(sprintf("  Interpretazione: %s\n", as.character(d_tg$magnitude)))
  cat("  (d positivo → femmine > maschi)\n")

  es_tg <- list(type = "Cohen's d", est = d_tg$estimate,
                ci = d_tg$conf.int, mag = as.character(d_tg$magnitude))
  test_tg <- list(name = "t-test", stat = t_tg$statistic,
                  df = t_tg$parameter, p = t_tg$p.value, ci = t_tg$conf.int)
} else {
  cat("→ Test non parametrico: Mann-Whitney U + Cliff's delta\n\n")
  w_tg <- wilcox.test(tg_f, tg_m, conf.int = TRUE)
  cat(sprintf("Mann-Whitney U: W = %.1f, p = %.4f\n",
              w_tg$statistic, w_tg$p.value))
  cat(sprintf("  Median maschi  = %.4f (IQR = %.4f, n = %d)\n",
              median(tg_m), IQR(tg_m), length(tg_m)))
  cat(sprintf("  Median femmine = %.4f (IQR = %.4f, n = %d)\n",
              median(tg_f), IQR(tg_f), length(tg_f)))
  cat(sprintf("  Hodges-Lehmann estimate = %.4f\n", w_tg$estimate))
  cat(sprintf("  IC 95%%: [%.4f, %.4f]\n",
              w_tg$conf.int[1], w_tg$conf.int[2]))
  cat(sprintf("\n  p < 0.025 (Bonferroni)? %s\n",
              ifelse(w_tg$p.value < 0.025, "SI", "NO")))

  cd_tg <- cliff.delta(tg_f, tg_m)
  cat(sprintf("\n  Cliff's delta = %.3f  [IC 95%%: %.3f, %.3f]\n",
              cd_tg$estimate, cd_tg$conf.int[1], cd_tg$conf.int[2]))
  cat(sprintf("  Interpretazione: %s\n", as.character(cd_tg$magnitude)))
  cat("  (delta positivo → femmine > maschi)\n")

  es_tg <- list(type = "Cliff's delta", est = cd_tg$estimate,
                ci = cd_tg$conf.int, mag = as.character(cd_tg$magnitude))
  test_tg <- list(name = "Mann-Whitney", stat = w_tg$statistic,
                  df = NA, p = w_tg$p.value, ci = w_tg$conf.int)
}

# --- CONFRONTO 2: WT maschi vs WT femmine ------------------------------------
cat("\n=======================================================================\n")
cat("  CONFRONTO 2: WT — MASCHI vs FEMMINE\n")
cat("  (Il dimorfismo sessuale è presente anche nei wild-type?)\n")
cat("=======================================================================\n\n")

wt_m <- long$value[long$genotype == "WT" & long$sex == "male"]
wt_f <- long$value[long$genotype == "WT" & long$sex == "female"]

norm_wt_m <- shapiro.test(wt_m)
norm_wt_f <- shapiro.test(wt_f)
wt_normal <- (norm_wt_m$p.value > 0.05) & (norm_wt_f$p.value > 0.05)

if (wt_normal) {
  cat("→ Test parametrico: Welch t-test + Cohen's d\n\n")
  t_wt <- t.test(wt_f, wt_m)
  cat(sprintf("Welch t-test: t = %.3f, df = %.1f, p = %.4f\n",
              t_wt$statistic, t_wt$parameter, t_wt$p.value))
  cat(sprintf("  Mean maschi  = %.4f (SD = %.4f, n = %d)\n",
              mean(wt_m), sd(wt_m), length(wt_m)))
  cat(sprintf("  Mean femmine = %.4f (SD = %.4f, n = %d)\n",
              mean(wt_f), sd(wt_f), length(wt_f)))
  cat(sprintf("  Differenza media (femmine - maschi) = %.4f\n",
              mean(wt_f) - mean(wt_m)))
  cat(sprintf("  IC 95%%: [%.4f, %.4f]\n",
              t_wt$conf.int[1], t_wt$conf.int[2]))
  cat(sprintf("\n  p < 0.025 (Bonferroni)? %s\n",
              ifelse(t_wt$p.value < 0.025, "SI", "NO")))

  d_wt <- cohen.d(wt_f, wt_m)
  cat(sprintf("\n  Cohen's d = %.3f  [IC 95%%: %.3f, %.3f]\n",
              d_wt$estimate, d_wt$conf.int[1], d_wt$conf.int[2]))
  cat(sprintf("  Interpretazione: %s\n", as.character(d_wt$magnitude)))
  cat("  (d positivo → femmine > maschi)\n")

  es_wt <- list(type = "Cohen's d", est = d_wt$estimate,
                ci = d_wt$conf.int, mag = as.character(d_wt$magnitude))
  test_wt <- list(name = "t-test", stat = t_wt$statistic,
                  df = t_wt$parameter, p = t_wt$p.value, ci = t_wt$conf.int)
} else {
  cat("→ Test non parametrico: Mann-Whitney U + Cliff's delta\n\n")
  w_wt <- wilcox.test(wt_f, wt_m, conf.int = TRUE)
  cat(sprintf("Mann-Whitney U: W = %.1f, p = %.4f\n",
              w_wt$statistic, w_wt$p.value))
  cat(sprintf("  Median maschi  = %.4f (IQR = %.4f, n = %d)\n",
              median(wt_m), IQR(wt_m), length(wt_m)))
  cat(sprintf("  Median femmine = %.4f (IQR = %.4f, n = %d)\n",
              median(wt_f), IQR(wt_f), length(wt_f)))
  cat(sprintf("  Hodges-Lehmann estimate = %.4f\n", w_wt$estimate))
  cat(sprintf("  IC 95%%: [%.4f, %.4f]\n",
              w_wt$conf.int[1], w_wt$conf.int[2]))
  cat(sprintf("\n  p < 0.025 (Bonferroni)? %s\n",
              ifelse(w_wt$p.value < 0.025, "SI", "NO")))

  cd_wt <- cliff.delta(wt_f, wt_m)
  cat(sprintf("\n  Cliff's delta = %.3f  [IC 95%%: %.3f, %.3f]\n",
              cd_wt$estimate, cd_wt$conf.int[1], cd_wt$conf.int[2]))
  cat(sprintf("  Interpretazione: %s\n", as.character(cd_wt$magnitude)))
  cat("  (delta positivo → femmine > maschi)\n")

  es_wt <- list(type = "Cliff's delta", est = cd_wt$estimate,
                ci = cd_wt$conf.int, mag = as.character(cd_wt$magnitude))
  test_wt <- list(name = "Mann-Whitney", stat = w_wt$statistic,
                  df = NA, p = w_wt$p.value, ci = w_wt$conf.int)
}

# --- Riepilogo ---------------------------------------------------------------
cat("\n=======================================================================\n")
cat("  RIEPILOGO — ATTIVITÀ ENZIMATICA FAAH\n")
cat("  Confronto: femmine vs maschi (vehicle)\n")
cat("=======================================================================\n\n")

cat(sprintf("%-20s %10s %8s %10s %10s %6s %s\n",
            "Confronto", "Test", "Stat", "p", "ES", "FC", "Bonf.(p<.025)"))
cat(strrep("-", 82), "\n")
cat(sprintf("%-20s %10s %8.3f %10.4f %10.3f %6.2f %s\n",
            "Tg (F vs M)", es_tg$type,
            test_tg$stat, test_tg$p, es_tg$est, fc_tg,
            ifelse(test_tg$p < 0.025, "SI *", "NO")))
cat(sprintf("%-20s %10s %8.3f %10.4f %10.3f %6.2f %s\n",
            "WT (F vs M)", es_wt$type,
            test_wt$stat, test_wt$p, es_wt$est, fc_wt,
            ifelse(test_wt$p < 0.025, "SI *", "NO")))

cat("\nFC = fold change (femmine / maschi). FC > 1 = attività maggiore nelle femmine.\n")
cat("ES positivo = femmine > maschi.\n")

# --- Chiudi sink -------------------------------------------------------------
sink()

# --- Salva oggetti -----------------------------------------------------------
save(long, desc, norm, model, test_tg, es_tg, test_wt, es_wt, fc_tg, fc_wt,
     file = "OUTPUT/08_FAAH_activity_sex_genotype.RData")

cat("\nFile salvati in OUTPUT/:\n")
cat(" -", basename(output_file), "\n")
cat(" - 08_FAAH_activity_sex_genotype.RData\n")
