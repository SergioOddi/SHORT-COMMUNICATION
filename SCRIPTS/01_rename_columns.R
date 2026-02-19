# =============================================================================
# 01_rename_columns.R
# Rinomina le colonne di tutti i dataset con convenzione inglese uniforme
# basata su Tortolani et al. 2025
# =============================================================================

library(readr)

# --- Dizionario di rinomina --------------------------------------------------
col_map <- c(
  "WT maschio veh"      = "WT_male_veh",
  "WT femmina veh"      = "WT_female_veh",
  "WT maschio PEA"      = "WT_male_PEA",
  "WT Femmina PEA"      = "WT_female_PEA",
  "Tg2576 maschio veh"  = "Tg_male_veh",
  "Tg2576 femmina veh"  = "Tg_female_veh",
  "Tg2576 maschio PEA"  = "Tg_male_PEA",
  "Tg2576 femmina PEA"  = "Tg_female_PEA"
)

# --- Funzione di rinomina ----------------------------------------------------
rename_and_save <- function(filename) {
  path <- file.path("DATASET", filename)
  df <- read_csv(path, show_col_types = FALSE)

  # Rinomina solo le colonne presenti nel dizionario
  old_names <- names(df)
  new_names <- ifelse(old_names %in% names(col_map), col_map[old_names], old_names)
  names(df) <- new_names

  write_csv(df, path)
  message(paste("Rinominato:", filename))
  message(paste("  Colonne:", paste(names(df), collapse = ", ")))
}

# --- Applica a tutti i dataset -----------------------------------------------
datasets <- c(
  "nor_test.csv",
  "apical_spine.csv",
  "basal_spine.csv",
  "3nt.csv",
  "iba1.csv",
  "gfap.csv",
  "inos.csv",
  "faah_genico.csv",
  "faah_proteica.csv",
  "naah_genico.csv"
)

invisible(lapply(datasets, rename_and_save))

message("\nRinomina completata per tutti i dataset.")
