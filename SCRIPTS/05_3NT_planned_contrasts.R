# =============================================================================
# 05_3NT_planned_contrasts.R
# Confronti pianificati su 3-Nitrotyrosine — Solo Tg2576
#
# Ipotesi a priori (basata su Tortolani et al., 2025):
#   La PEA riduce lo stress nitrosativo nei Tg2576.
#   Ci aspettiamo una riduzione in entrambi i sessi.
#
# Confronti:
#   1. Tg_male_veh   vs Tg_male_PEA   (replicazione dato pubblicato)
#   2. Tg_female_veh vs Tg_female_PEA  (estensione alle femmine)
#
# Tutti i gruppi normali (Shapiro-Wilk p > 0.05) → test parametrici
# Correzione: Bonferroni per 2 confronti (α = 0.025)
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
df <- read_csv("DATASET/3nt.csv", show_col_types = FALSE)

long <- df |>
  pivot_longer(cols = everything(), names_to = "group", values_to = "value") |>
  filter(!is.na(value)) |>
  mutate(
    sex       = ifelse(grepl("female", group), "female", "male"),
    treatment = ifelse(grepl("PEA", group), "PEA", "veh")
  )

# --- Output ------------------------------------------------------------------
if (!dir.exists("OUTPUT")) dir.create("OUTPUT")

output_file <- "OUTPUT/05_3NT_planned_contrasts.txt"
sink(output_file, split = TRUE)

cat("=======================================================================\n")
cat("  3-NITROTYROSINE — CONFRONTI PIANIFICATI (Solo Tg2576)\n")
cat("  Ipotesi a priori: PEA riduce lo stress nitrosativo in entrambi i sessi\n")
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

# --- Test di normalità -------------------------------------------------------
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

# --- ANOVA 2 vie di contesto ------------------------------------------------
cat("\n=======================================================================\n")
cat("  ANOVA 2 VIE (Sex × Treatment) — analisi omnibus di contesto\n")
cat("=======================================================================\n\n")

model <- aov(value ~ sex * treatment, data = long)
print(summary(model))

ss <- summary(model)[[1]]
ss_total <- sum(ss[["Sum Sq"]])
cat("\nEta-squared (η²):\n")
cat(sprintf("  Sex:              η² = %.4f\n", ss[["Sum Sq"]][1] / ss_total))
cat(sprintf("  Treatment:        η² = %.4f\n", ss[["Sum Sq"]][2] / ss_total))
cat(sprintf("  Sex × Treatment:  η² = %.4f\n", ss[["Sum Sq"]][3] / ss_total))

# --- CONFRONTO 1: Maschi (veh vs PEA) ---------------------------------------
cat("\n=======================================================================\n")
cat("  CONFRONTO 1: Tg MASCHI — veh vs PEA\n")
cat("  (Replicazione del dato pubblicato, Tortolani et al. 2025)\n")
cat("  Nota: valori più bassi = meno stress nitrosativo = effetto positivo\n")
cat("=======================================================================\n\n")

m_veh <- long$value[long$sex == "male" & long$treatment == "veh"]
m_pea <- long$value[long$sex == "male" & long$treatment == "PEA"]

t_male <- t.test(m_veh, m_pea)
cat(sprintf("Welch t-test: t = %.3f, df = %.1f, p = %.4f\n",
            t_male$statistic, t_male$parameter, t_male$p.value))
cat(sprintf("  Mean veh = %.2f (SD = %.2f, n = %d)\n",
            mean(m_veh), sd(m_veh), length(m_veh)))
cat(sprintf("  Mean PEA = %.2f (SD = %.2f, n = %d)\n",
            mean(m_pea), sd(m_pea), length(m_pea)))
cat(sprintf("  Differenza media (veh - PEA) = %.2f\n", mean(m_veh) - mean(m_pea)))
cat(sprintf("  IC 95%% della differenza: [%.2f, %.2f]\n",
            t_male$conf.int[1], t_male$conf.int[2]))

cat(sprintf("\n  p < 0.025 (Bonferroni)? %s\n",
            ifelse(t_male$p.value < 0.025, "SI", "NO")))

# Cohen's d (veh vs PEA: d positivo = veh > PEA = PEA riduce 3-NT)
d_male <- cohen.d(m_veh, m_pea)
cat(sprintf("\n  Cohen's d = %.3f  [IC 95%%: %.3f, %.3f]\n",
            d_male$estimate, d_male$conf.int[1], d_male$conf.int[2]))
cat(sprintf("  Interpretazione: %s\n", as.character(d_male$magnitude)))
cat("  (d positivo → veh > PEA → PEA riduce 3-NT)\n")

# --- CONFRONTO 2: Femmine (veh vs PEA) --------------------------------------
cat("\n=======================================================================\n")
cat("  CONFRONTO 2: Tg FEMMINE — veh vs PEA\n")
cat("  (Estensione alle femmine — dato nuovo)\n")
cat("  Nota: valori più bassi = meno stress nitrosativo = effetto positivo\n")
cat("=======================================================================\n\n")

f_veh <- long$value[long$sex == "female" & long$treatment == "veh"]
f_pea <- long$value[long$sex == "female" & long$treatment == "PEA"]

t_female <- t.test(f_veh, f_pea)
cat(sprintf("Welch t-test: t = %.3f, df = %.1f, p = %.4f\n",
            t_female$statistic, t_female$parameter, t_female$p.value))
cat(sprintf("  Mean veh = %.2f (SD = %.2f, n = %d)\n",
            mean(f_veh), sd(f_veh), length(f_veh)))
cat(sprintf("  Mean PEA = %.2f (SD = %.2f, n = %d)\n",
            mean(f_pea), sd(f_pea), length(f_pea)))
cat(sprintf("  Differenza media (veh - PEA) = %.2f\n", mean(f_veh) - mean(f_pea)))
cat(sprintf("  IC 95%% della differenza: [%.2f, %.2f]\n",
            t_female$conf.int[1], t_female$conf.int[2]))

cat(sprintf("\n  p < 0.025 (Bonferroni)? %s\n",
            ifelse(t_female$p.value < 0.025, "SI", "NO")))

d_female <- cohen.d(f_veh, f_pea)
cat(sprintf("\n  Cohen's d = %.3f  [IC 95%%: %.3f, %.3f]\n",
            d_female$estimate, d_female$conf.int[1], d_female$conf.int[2]))
cat(sprintf("  Interpretazione: %s\n", as.character(d_female$magnitude)))
cat("  (d positivo → veh > PEA → PEA riduce 3-NT)\n")

# --- Riepilogo ---------------------------------------------------------------
cat("\n=======================================================================\n")
cat("  RIEPILOGO — 3-NITROTYROSINE\n")
cat("=======================================================================\n\n")

cat(sprintf("%-20s %8s %8s %10s %10s %s\n",
            "Confronto", "t", "p", "Cohen's d", "Magnitude", "Bonf.(p<.025)"))
cat(strrep("-", 75), "\n")
cat(sprintf("%-20s %8.3f %8.4f %10.3f %10s %s\n",
            "Maschi (veh vs PEA)",
            t_male$statistic, t_male$p.value,
            d_male$estimate, as.character(d_male$magnitude),
            ifelse(t_male$p.value < 0.025, "SI *", "NO")))
cat(sprintf("%-20s %8.3f %8.4f %10.3f %10s %s\n",
            "Femmine (veh vs PEA)",
            t_female$statistic, t_female$p.value,
            d_female$estimate, as.character(d_female$magnitude),
            ifelse(t_female$p.value < 0.025, "SI *", "NO")))

cat("\n-----------------------------------------------------------------------\n")
cat("Interpretazione Cohen's d:\n")
cat("  |d| < 0.2 = trascurabile, 0.2-0.5 = piccolo,\n")
cat("  0.5-0.8 = medio, > 0.8 = grande\n")
cat("  d positivo = veh > PEA (PEA riduce stress nitrosativo)\n")
cat("-----------------------------------------------------------------------\n")

# --- Chiudi sink -------------------------------------------------------------
sink()

# --- Salva oggetti -----------------------------------------------------------
save(long, desc, norm, model, t_male, t_female, d_male, d_female,
     file = "OUTPUT/05_3NT_planned_contrasts.RData")

cat("\nFile salvati in OUTPUT/:\n")
cat(" -", basename(output_file), "\n")
cat(" - 05_3NT_planned_contrasts.RData\n")
