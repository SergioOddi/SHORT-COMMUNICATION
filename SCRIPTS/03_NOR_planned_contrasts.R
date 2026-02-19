# =============================================================================
# 03_NOR_planned_contrasts.R
# Confronti pianificati sul NOR Test — Solo Tg2576
#
# Ipotesi a priori (basata su Tortolani et al., 2025):
#   La PEA migliora la performance cognitiva nei Tg2576.
#   Ci aspettiamo un effetto positivo sia nei maschi che nelle femmine.
#
# Confronti:
#   1. Tg_male_veh   vs Tg_male_PEA   (replicazione dato pubblicato)
#   2. Tg_female_veh vs Tg_female_PEA  (estensione alle femmine)
#
# Correzione: Bonferroni per 2 confronti (α = 0.025 per confronto)
# Effect size: Cohen's d con IC 95%
# =============================================================================

# --- Pacchetti ---------------------------------------------------------------
required_pkgs <- c("readr", "tidyr", "dplyr", "effsize")
new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs)) install.packages(new_pkgs)

library(readr)
library(tidyr)
library(dplyr)
library(effsize)

# --- Caricamento dati --------------------------------------------------------
df <- read_csv("DATASET/nor_test.csv", show_col_types = FALSE)

# Converti in formato long e filtra solo Tg2576
long <- df |>
  pivot_longer(cols = everything(), names_to = "group", values_to = "value") |>
  filter(!is.na(value)) |>
  filter(grepl("^Tg", group)) |>
  mutate(
    sex       = ifelse(grepl("female", group), "female", "male"),
    treatment = ifelse(grepl("PEA", group), "PEA", "veh")
  )

# --- Crea cartella OUTPUT se non esiste --------------------------------------
if (!dir.exists("OUTPUT")) dir.create("OUTPUT")

# --- Avvia cattura output ----------------------------------------------------
output_file <- "OUTPUT/03_NOR_planned_contrasts.txt"
sink(output_file, split = TRUE)

cat("=======================================================================\n")
cat("  NOR TEST — CONFRONTI PIANIFICATI (Solo Tg2576)\n")
cat("  Ipotesi a priori: PEA migliora la performance in entrambi i sessi\n")
cat("  Correzione Bonferroni per 2 confronti (α = 0.025)\n")
cat("=======================================================================\n")
cat("Analisi eseguita il:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# --- Statistiche descrittive -------------------------------------------------
cat("\n--- Statistiche descrittive ---\n\n")
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

# --- Test di normalità per sottogruppo ---------------------------------------
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
cat("\nTutti i gruppi normali (p > 0.05): ", all(norm$SW_p > 0.05), "\n")

# --- ANOVA 2 vie di contesto (Sex × Treatment) ------------------------------
cat("\n=======================================================================\n")
cat("  ANOVA 2 VIE (Sex × Treatment) — analisi omnibus di contesto\n")
cat("=======================================================================\n\n")

model <- aov(value ~ sex * treatment, data = long)
print(summary(model))

# Eta-squared (effect size per ANOVA)
ss <- summary(model)[[1]]
ss_total <- sum(ss[["Sum Sq"]])
cat("\nEta-squared (η²):\n")
cat(sprintf("  Sex:              η² = %.4f\n", ss[["Sum Sq"]][1] / ss_total))
cat(sprintf("  Treatment:        η² = %.4f\n", ss[["Sum Sq"]][2] / ss_total))
cat(sprintf("  Sex × Treatment:  η² = %.4f\n", ss[["Sum Sq"]][3] / ss_total))

# --- CONFRONTO 1: Maschi Tg (veh vs PEA) ------------------------------------
cat("\n=======================================================================\n")
cat("  CONFRONTO 1: Tg MASCHI — veh vs PEA\n")
cat("  (Replicazione del dato pubblicato, Tortolani et al. 2025)\n")
cat("=======================================================================\n\n")

m_veh <- long$value[long$sex == "male" & long$treatment == "veh"]
m_pea <- long$value[long$sex == "male" & long$treatment == "PEA"]

# Welch t-test
t_male <- t.test(m_pea, m_veh)
cat(sprintf("Welch t-test: t = %.3f, df = %.1f, p = %.4f\n",
            t_male$statistic, t_male$parameter, t_male$p.value))
cat(sprintf("  Mean veh = %.2f (SD = %.2f, n = %d)\n",
            mean(m_veh), sd(m_veh), length(m_veh)))
cat(sprintf("  Mean PEA = %.2f (SD = %.2f, n = %d)\n",
            mean(m_pea), sd(m_pea), length(m_pea)))
cat(sprintf("  Differenza media = %.2f\n", mean(m_pea) - mean(m_veh)))
cat(sprintf("  IC 95%% della differenza: [%.2f, %.2f]\n",
            t_male$conf.int[1], t_male$conf.int[2]))

# Significativo dopo Bonferroni?
cat(sprintf("\n  p < 0.025 (Bonferroni)? %s\n",
            ifelse(t_male$p.value < 0.025, "SI", "NO")))

# Cohen's d
d_male <- cohen.d(m_pea, m_veh)
cat(sprintf("\n  Cohen's d = %.3f  [IC 95%%: %.3f, %.3f]\n",
            d_male$estimate, d_male$conf.int[1], d_male$conf.int[2]))
cat(sprintf("  Interpretazione: %s\n", as.character(d_male$magnitude)))

# --- CONFRONTO 2: Femmine Tg (veh vs PEA) -----------------------------------
cat("\n=======================================================================\n")
cat("  CONFRONTO 2: Tg FEMMINE — veh vs PEA\n")
cat("  (Estensione alle femmine — dato nuovo)\n")
cat("=======================================================================\n\n")

f_veh <- long$value[long$sex == "female" & long$treatment == "veh"]
f_pea <- long$value[long$sex == "female" & long$treatment == "PEA"]

# Welch t-test
t_female <- t.test(f_pea, f_veh)
cat(sprintf("Welch t-test: t = %.3f, df = %.1f, p = %.4f\n",
            t_female$statistic, t_female$parameter, t_female$p.value))
cat(sprintf("  Mean veh = %.2f (SD = %.2f, n = %d)\n",
            mean(f_veh), sd(f_veh), length(f_veh)))
cat(sprintf("  Mean PEA = %.2f (SD = %.2f, n = %d)\n",
            mean(f_pea), sd(f_pea), length(f_pea)))
cat(sprintf("  Differenza media = %.2f\n", mean(f_pea) - mean(f_veh)))
cat(sprintf("  IC 95%% della differenza: [%.2f, %.2f]\n",
            t_female$conf.int[1], t_female$conf.int[2]))

# Significativo dopo Bonferroni?
cat(sprintf("\n  p < 0.025 (Bonferroni)? %s\n",
            ifelse(t_female$p.value < 0.025, "SI", "NO")))

# Cohen's d
d_female <- cohen.d(f_pea, f_veh)
cat(sprintf("\n  Cohen's d = %.3f  [IC 95%%: %.3f, %.3f]\n",
            d_female$estimate, d_female$conf.int[1], d_female$conf.int[2]))
cat(sprintf("  Interpretazione: %s\n", as.character(d_female$magnitude)))

# --- Riepilogo confronti ----------------------------------------------------
cat("\n=======================================================================\n")
cat("  RIEPILOGO\n")
cat("=======================================================================\n\n")

cat(sprintf("%-20s %8s %8s %10s %10s %s\n",
            "Confronto", "t", "p", "Cohen's d", "Magnitude", "Bonf.(p<.025)"))
cat(strrep("-", 75), "\n")
cat(sprintf("%-20s %8.3f %8.4f %10.3f %10s %s\n",
            "Maschi (PEA vs veh)",
            t_male$statistic, t_male$p.value,
            d_male$estimate, as.character(d_male$magnitude),
            ifelse(t_male$p.value < 0.025, "SI *", "NO")))
cat(sprintf("%-20s %8.3f %8.4f %10.3f %10s %s\n",
            "Femmine (PEA vs veh)",
            t_female$statistic, t_female$p.value,
            d_female$estimate, as.character(d_female$magnitude),
            ifelse(t_female$p.value < 0.025, "SI *", "NO")))

cat("\n-----------------------------------------------------------------------\n")
cat("Interpretazione:\n")
cat("  Cohen's d: |d| < 0.2 = trascurabile, 0.2-0.5 = piccolo,\n")
cat("             0.5-0.8 = medio, > 0.8 = grande\n")
cat("-----------------------------------------------------------------------\n")

# --- Chiudi cattura output ---------------------------------------------------
sink()

# --- Salva oggetti -----------------------------------------------------------
save(long, desc, norm, model, t_male, t_female, d_male, d_female,
     file = "OUTPUT/03_NOR_planned_contrasts.RData")

cat("\nFile salvati in OUTPUT/:\n")
cat(" -", basename(output_file), "\n")
cat(" - 03_NOR_planned_contrasts.RData\n")
