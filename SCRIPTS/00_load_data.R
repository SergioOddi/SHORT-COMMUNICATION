# =============================================================================
# 00_load_data.R
# Scarica il dataset da Google Sheets e lo salva in DATASET/ come CSV
# =============================================================================

# --- Pacchetti richiesti -----------------------------------------------------
if (!requireNamespace("googlesheets4", quietly = TRUE)) {
  install.packages("googlesheets4")
}
if (!requireNamespace("readr", quietly = TRUE)) {
  install.packages("readr")
}

library(googlesheets4)
library(readr)

# --- Autenticazione ----------------------------------------------------------
# La prima volta apre il browser per il login Google.
# Le credenziali vengono memorizzate in cache per le sessioni successive.
gs4_auth()

# --- URL del foglio Google Sheets --------------------------------------------
# Sostituisci con l'URL del tuo foglio (o solo l'ID del foglio)
SHEET_URL <- "INSERISCI_QUI_URL_DEL_FOGLIO_GOOGLE_SHEETS"

# Opzionale: specifica il nome o l'indice del foglio (tab) se ce ne sono più
# SHEET_NAME <- "Foglio1"  # decommenta se necessario

# --- Caricamento dati --------------------------------------------------------
message("Connessione a Google Sheets in corso...")

df_raw <- read_sheet(SHEET_URL)
# Con foglio specifico: df_raw <- read_sheet(SHEET_URL, sheet = SHEET_NAME)

message(paste("Dataset caricato:", nrow(df_raw), "righe x", ncol(df_raw), "colonne"))
message("Colonne: ", paste(names(df_raw), collapse = ", "))

# --- Salvataggio in CSV ------------------------------------------------------
output_path <- "DATASET/dataset_raw.csv"

write_csv(df_raw, output_path)

message(paste("Dataset salvato in:", output_path))

# --- Anteprima ---------------------------------------------------------------
print(head(df_raw))
