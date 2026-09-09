rm(list = ls())
library(readxl)

mainDir <- "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-case-study-mm"
input_file <- file.path(mainDir, "data/EloKRd cfDNA PET raw data_6_2025.xlsx")

# ============================================================
# 1. READ RAW DATA
# ============================================================
raw <- read_excel(input_file, col_names = FALSE)

# Row 1 = Excel metadata, Row 2 = column names
real_col_names <- trimws(as.character(raw[2, ]))
raw_data <- raw[-(1:2), ]
colnames(raw_data) <- real_col_names
raw_data <- as.data.frame(raw_data, stringsAsFactors = FALSE, check.names = FALSE)

# Drop trailing blank rows
id_col <- which(colnames(raw_data) == "Subject ID")
if (length(id_col) > 0) {
  raw_data <- raw_data[!is.na(raw_data[[id_col]]) & raw_data[[id_col]] != "", ]
}
cat(sprintf("EloKRD raw: %d subjects x %d columns\n", nrow(raw_data), ncol(raw_data)))

# ============================================================
# 2. HANDLE DUPLICATE COLUMNS (keep the one with fewest NAs)
# ============================================================
col_names <- colnames(raw_data)
dup_names <- unique(col_names[duplicated(col_names)])
cols_to_drop <- c()

for (dn in dup_names) {
  idx <- which(col_names == dn)
  na_counts <- sapply(idx, function(i) sum(is.na(raw_data[[i]]) | raw_data[[i]] == ""))
  min_na <- min(na_counts)
  candidates <- idx[na_counts == min_na]
  keep <- candidates[length(candidates)]
  drop <- setdiff(idx, keep)
  cat(sprintf('  Duplicate "%s": keeping col %d, dropping col %s\n',
              dn, keep, paste(drop, collapse = ", ")))
  cols_to_drop <- c(cols_to_drop, drop)
}

# Drop Notes column if in expected position
notes_col <- which(col_names == "Notes")
if (length(notes_col) == 1) {
  cols_to_drop <- c(cols_to_drop, notes_col)
}

if (length(cols_to_drop) > 0) {
  raw_data <- raw_data[, -cols_to_drop]
}
cat(sprintf("After dedup: %d columns\n", ncol(raw_data)))

# ============================================================
# 3. EXTRACT AND HARMONIZE ANALYSIS VARIABLES
# ============================================================
# Helper: safely get column (returns NA vector if missing)
# Uses fuzzy matching to handle extra whitespace in column names
get_col <- function(df, name) {
  if (name %in% names(df)) return(df[[name]])
  # Try collapsing whitespace for matching
  norm <- function(x) gsub("\\s+", " ", trimws(x))
  idx <- which(norm(names(df)) == norm(name))
  if (length(idx) == 1) return(df[[idx]])
  rep(NA, nrow(df))
}

to_numeric <- function(x) suppressWarnings(as.numeric(x))

elokrd <- data.frame(
  # ID (synthetic, sequential -- NOT the institutional Subject ID)
  subject_id     = sprintf("EloKRd-%02d", seq_len(nrow(raw_data))),
  source         = "EloKRD",

  # Demographics
  age            = to_numeric(get_col(raw_data, "Age at Enrollment")),
  sex            = get_col(raw_data, "Sex"),
  race           = get_col(raw_data, "Race"),
  ethnicity      = get_col(raw_data, "Ethnicity"),

  # Disease staging
  iss            = get_col(raw_data, "ISS staging"),
  riss           = get_col(raw_data, "R-ISS staging"),
  ecog           = get_col(raw_data, "ECOG at time of enrollment"),

  # Cytogenetics / FISH
  high_risk_cyto = to_numeric(get_col(raw_data, "High-risk Cyto (1=Yes, 0 = No)")),
  t_4_14         = to_numeric(get_col(raw_data, "t(4, 14) abnormality")),
  t_14_16        = to_numeric(get_col(raw_data, "t(14, 16) abnormality")),
  t_14_20        = to_numeric(get_col(raw_data, "t(14, 20) abnormality")),
  del_17p        = to_numeric(get_col(raw_data, "17p13 del abnormality")),
  gain_1q        = to_numeric(get_col(raw_data, "1q21 gain abnormality")),
  del_1p         = to_numeric(get_col(raw_data, "del 1p abnormality")),
  double_hit     = get_col(raw_data, "Double-hit myeloma?"),

  # Heavy/light chains
  heavy_chain    = get_col(raw_data, "Heavy Chains (IgG or IgA)"),
  light_chain    = get_col(raw_data, "Light Chains (Kappa or Lambda)"),

  # Treatment details
  transplant_off_protocol = to_numeric(get_col(raw_data,
                           "Transplant off-protocol (1 = Yes, 0 = No)")),
  best_response  = get_col(raw_data, "Best Overall Response through Cutoff"),

  # Survival outcomes (from C1D1)
  os_months      = to_numeric(get_col(raw_data, "OS through Cut-off Date (months)")),
  os_status      = to_numeric(get_col(raw_data, "OS status")),
  pfs_months     = to_numeric(get_col(raw_data, "PFS (months)")),
  # Source uses 2 = death without prior progression; recode to 1 (still a PFS event)
  pfs_status     = ifelse(to_numeric(get_col(raw_data, "PFS status (1 = prog or death)")) == 2, 1,
                          to_numeric(get_col(raw_data, "PFS status (1 = prog or death)"))),

  # NOTE: calendar dates (C1D1/enrollment/dx/progression/death/last follow-up)
  # are intentionally NOT retained here -- only de-identified relative times
  # (os_months, pfs_months) are kept for public release.

  stringsAsFactors = FALSE
)

# ============================================================
# 4. RESTRICT TO NO-ASCT PATIENTS (dropped — ASCT kept as covariate)
# ============================================================
# n_before <- nrow(elokrd)
# elokrd <- elokrd[elokrd$transplant_off_protocol == 0 & !is.na(elokrd$transplant_off_protocol), ]
# cat(sprintf("\nRestricted to no-ASCT: %d -> %d (dropped %d with ASCT)\n",
#             n_before, nrow(elokrd), n_before - nrow(elokrd)))
cat(sprintf("\nASCT (transplant_off_protocol) retained as covariate: %d with ASCT, %d without\n",
            sum(elokrd$transplant_off_protocol == 1, na.rm = TRUE),
            sum(elokrd$transplant_off_protocol == 0, na.rm = TRUE)))

# ============================================================
# 5. VALIDATE
# ============================================================
cat(sprintf("\n=== EloKRD Cleaned: %d subjects x %d columns ===\n",
            nrow(elokrd), ncol(elokrd)))
cat(sprintf("  OS events:  %d / %d\n", sum(elokrd$os_status == 1, na.rm = TRUE), nrow(elokrd)))
cat(sprintf("  PFS events: %d / %d\n", sum(elokrd$pfs_status == 1, na.rm = TRUE), nrow(elokrd)))
cat(sprintf("  Median OS:  %.1f months\n", median(elokrd$os_months, na.rm = TRUE)))
cat(sprintf("  Median PFS: %.1f months\n", median(elokrd$pfs_months, na.rm = TRUE)))

cat("\n  Sex distribution:\n")
print(table(elokrd$sex, useNA = "ifany"))
cat("\n  Race distribution:\n")
print(table(elokrd$race, useNA = "ifany"))
cat("\n  ISS distribution:\n")
print(table(elokrd$iss, useNA = "ifany"))
cat("\n  High-risk cytogenetics:\n")
print(table(elokrd$high_risk_cyto, useNA = "ifany"))

# ============================================================
# 6. SAVE TO data_cleaned/
# ============================================================
out_dir <- file.path(mainDir, "data_cleaned")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_file <- file.path(out_dir, sprintf("elokrd_cleaned_n%d.RData", nrow(elokrd)))
save(elokrd, file = out_file)
cat(sprintf("\nSaved: %s\n", out_file))
