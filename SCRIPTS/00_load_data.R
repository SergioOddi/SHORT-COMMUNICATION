# =============================================================================
# 00_load_data.R
# Scarica i 10 sheet dal file Google Sheets e li salva in DATASET/ come CSV
# =============================================================================

# --- Pacchetti richiesti -----------------------------------------------------
if (!requireNamespace("googlesheets4", quietly = TRUE)) install.packages("googlesheets4")
if (!requireNamespace("readr", quietly = TRUE))        install.packages("readr")

library(googlesheets4)
library(readr)

# --- Autenticazione ----------------------------------------------------------
# La prima volta apre il browser per il login Google.
# Le credenziali vengono memorizzate in cache per le sessioni successive.
gs4_auth(scopes = "https://www.googleapis.com/auth/spreadsheets.readonly")

# --- URL del file Google Sheets ----------------------------------------------
SHEET_URL <- "https://docs.google.com/spreadsheets/d/1lxIAOL4Mv7LVprqps2IRHu9s6iK3hp4ajMOoYuexQrQ/edit?gid=0#gid=0"

# --- Nomi dei 10 sheet (tab) -------------------------------------------------
sheet_names <- list(
  nor_test      = "NOR Test",
  apical_spine  = "Apical Spine density",
  basal_spine   = "Basal Spine density",
  nt3           = "3-NT",
  iba1          = "IBA1",
  gfap          = "GFAP",
  inos          = "iNOS",
  faah_genico   = "FAAH genico",
  faah_proteica = "FAAH proteica",
  naah_genico   = "NAAH genico"
)

# --- Caricamento dei dataset -------------------------------------------------
load_sheet <- function(sheet_id, sheet_name) {
  message(paste("Caricamento:", sheet_name, "..."))
  df <- read_sheet(SHEET_URL, sheet = sheet_name)
  message(paste("  ->", nrow(df), "righe x", ncol(df), "colonne"))
  return(df)
}

datasets <- mapply(
  load_sheet,
  sheet_id   = names(sheet_names),
  sheet_name = unlist(sheet_names),
  SIMPLIFY   = FALSE
)

# Assegna ogni dataset a una variabile con nome leggibile
nor_test      <- datasets$nor_test
apical_spine  <- datasets$apical_spine
basal_spine   <- datasets$basal_spine
nt3           <- datasets$nt3
iba1          <- datasets$iba1
gfap          <- datasets$gfap
inos          <- datasets$inos
faah_genico   <- datasets$faah_genico
faah_proteica <- datasets$faah_proteica
naah_genico   <- datasets$naah_genico

# --- Salvataggio in CSV ------------------------------------------------------
message("\nSalvataggio dei CSV in DATASET/ ...")

save_csv <- function(df, filename) {
  path <- file.path("DATASET", filename)
  write_csv(df, path)
  message(paste("  Salvato:", path))
}

save_csv(nor_test,      "nor_test.csv")
save_csv(apical_spine,  "apical_spine.csv")
save_csv(basal_spine,   "basal_spine.csv")
save_csv(nt3,           "3nt.csv")
save_csv(iba1,          "iba1.csv")
save_csv(gfap,          "gfap.csv")
save_csv(inos,          "inos.csv")
save_csv(faah_genico,   "faah_genico.csv")
save_csv(faah_proteica, "faah_proteica.csv")
save_csv(naah_genico,   "naah_genico.csv")

message("\nTutti i dataset sono stati scaricati e salvati in DATASET/")
