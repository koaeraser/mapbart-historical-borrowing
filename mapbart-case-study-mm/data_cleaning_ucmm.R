rm(list = ls())

mainDir <- "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-case-study-mm"

# ============================================================
# 1. LOAD UCMM FULL DATA (n=797)
# ============================================================
load(file.path(mainDir, "data/UCMM_cleaned_full_n797.RData"))
ucmm <- analysis_data
rm(analysis_data)
cat(sprintf("UCMM full: %d subjects x %d columns\n", nrow(ucmm), ncol(ucmm)))

to_days   <- function(x) as.numeric(as.Date(x))
to_months <- function(days) days / 30.44

# ============================================================
# 2. COHORT RESTRICTIONS (applied IN ORDER; resulting n recorded below)
# ------------------------------------------------------------
#    Full UCMM (analysis_data) ............................ n = 797
# ------------------------------------------------------------
# Restriction order: 1) negative OS/PFS, 2) Has Tx, 3) ASCT on/before
# induction, 4) E-Rd.  OS/PFS are measured from the induction (Tx) start date
# (analogous to the EloKRd C1D1 date); they are computed on the full file first
# (NA when no Tx date) so the negative-time filter precedes the Has-Tx filter.
# ============================================================
# t=0 induction/Tx start; 
# NOT documented
tx_start_all <- to_days(ucmm$induction_time_abstracted)
# Death; 
# codebook: "Date of death from EMR"
dod_all      <- to_days(ucmm$dod_emr)
# Last day of contact; 
# codebook: Derived = 1) dod_emr, 2) later of last_contact_cr & last_encounter_emr
last_all     <- to_days(ucmm$last_contact_composite)
# PFS primary progression; 
# NOT documented
prog_all_d   <- to_days(ucmm$first_progression_feb2022_abstracted)
# PFS fallback; 
# NOT documented
recur_all_d  <- to_days(ucmm$recurrence_date_cr)

# Use first_progression_feb2022_abstracted as primary progression date;
# fall back to recurrence_date_cr if the former is missing.
prog_all <- ifelse(!is.na(prog_all_d), prog_all_d, recur_all_d)  

# OS event = death observed; OS time = (death | last contact) - induction
os_event_all  <- ifelse(!is.na(dod_all), 1, 0)
os_months_all <- ifelse(!is.na(dod_all), to_months(dod_all - tx_start_all),
                                          to_months(last_all - tx_start_all))
# PFS event = progression OR death; PFS time = first of (progression, death),
# else last contact -- all measured from induction (tx_start_all)
pfs_event_all <- ifelse(!is.na(prog_all) | !is.na(dod_all), 1, 0)
pfs_raw_all   <- pmin(ifelse(!is.na(prog_all), prog_all - tx_start_all, Inf),
                      ifelse(!is.na(dod_all),  dod_all  - tx_start_all, Inf), na.rm = TRUE)
pfs_months_all <- ifelse(is.finite(pfs_raw_all), to_months(pfs_raw_all),
                                                  to_months(last_all - tx_start_all))

# Restriction masks (TRUE = keep), applied cumulatively in the order above:
keep_neg  <- !((os_months_all < 0) %in% TRUE | (pfs_months_all < 0) %in% TRUE)  # Step 1
has_tx    <- !is.na(tx_start_all)                                               # Step 2
asct_rel  <- to_days(ucmm$asct_time_abstracted) - tx_start_all
keep_asct <- !((!is.na(asct_rel)) & asct_rel <= 0)                             # Step 3
keep_erd  <- !(as.character(ucmm$regimen) %in% c("E-Rd", "Elo-Rd"))            # Step 4

m1 <- keep_neg
m2 <- m1 & has_tx
m3 <- m2 & keep_asct
m4 <- m3 & keep_erd
cat(sprintf("Full UCMM: n=%d\n", nrow(ucmm)))
cat(sprintf("1) Exclude negative OS/PFS time (dropped %d): n=%d\n",    sum(!m1), sum(m1)))
cat(sprintf("2) Has Tx (dropped %d): n=%d\n",                         sum(m1 & !has_tx), sum(m2)))
cat(sprintf("3) Exclude ASCT on/before induction (dropped %d): n=%d\n",sum(m2 & !keep_asct), sum(m3)))
cat(sprintf("4) Exclude E-Rd (dropped %d): n=%d\n",                   sum(m3 & !keep_erd), sum(m4)))
#    (1) Exclude negative OS/PFS time (dropped 7) ......... n = 790
#    (2) Has Tx (dropped 487) ............................. n = 303   <- COHORT 1 (Has Tx)
#    (3) Exclude ASCT on/before induction (dropped 49) .... n = 254   <- COHORT 2 (excl. ASCT on/before)
#    (4) Exclude E-Rd (dropped 1) ......................... n = 253   <- PRIMARY analysis cohort (all regimens)

# PRIMARY cohort = all patients passing steps 1-4 (any regimen; E-Rd excluded)
idx        <- which(m4)
ucmm_sub   <- ucmm[idx, ]
os_event   <- os_event_all[idx];  os_months  <- os_months_all[idx]
pfs_event  <- pfs_event_all[idx]; pfs_months <- pfs_months_all[idx]
# Floor zero PFS times (same-day progression) at 0.1 month so log() is finite downstream
pfs_months <- pmax(pfs_months, 0.1)

# ============================================================
# 3. HELPER: safely get column
# ============================================================
get_col <- function(df, name) {
  if (name %in% names(df)) df[[name]] else rep(NA, nrow(df))
}

to_numeric_safe <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

# ============================================================
# 4. EXTRACT AND HARMONIZE ANALYSIS VARIABLES (on the primary cohort)
# ============================================================

# --- Sex ---
# Harmonize sex_errc values
sex_raw <- as.character(ucmm_sub$sex_errc)
sex_clean <- ifelse(grepl("[Ff]emale|^F$|^2$", sex_raw), "Female",
             ifelse(grepl("[Mm]ale|^M$|^1$", sex_raw), "Male", NA))

# --- Race ---
# Use race_composite (the harmonized version)
race_raw <- as.character(get_col(ucmm_sub, "race_composite"))

# --- Age ---
# Use age at identification date: compute from dob_emr and date_identified_errc
dob <- as.Date(ucmm_sub$dob_emr)
tx_date <- as.Date(ucmm_sub$date_identified_errc)
age_at_tx <- as.numeric(difftime(tx_date, dob, units = "days")) / 365.25
# Fall back to age_dx_cr if dob-based computation fails
age_dx <- to_numeric_safe(get_col(ucmm_sub, "age_dx_cr"))
age_final <- ifelse(!is.na(age_at_tx), age_at_tx, age_dx)

# --- ISS ---
# Use iss_derived (harmonized), fall back to iss_cr
iss_derived <- as.character(get_col(ucmm_sub, "iss_derived"))
iss_cr      <- as.character(get_col(ucmm_sub, "iss_cr"))
iss_final   <- ifelse(!is.na(iss_derived) & iss_derived != "", iss_derived, iss_cr)

# --- R-ISS ---
riss_derived <- as.character(get_col(ucmm_sub, "riss_derived"))

# --- Cytogenetics / FISH (use "at dx" variants, closest to baseline) ---
# These should be 0/1 or "Positive"/"Negative" etc. Harmonize to 0/1.
harmonize_cyto <- function(x) {
  x <- tolower(as.character(x))
  ifelse(x %in% c("1", "positive", "yes", "true", "abnormal"), 1,
  ifelse(x %in% c("0", "negative", "no", "false", "normal", "not detected"), 0, NA))
}

t_4_14  <- harmonize_cyto(get_col(ucmm_sub, "hierarchy_corrected_t(4:14) - FGFR3/IGH fusion at dx"))
# Fall back to non-hierarchy version
if (all(is.na(t_4_14))) t_4_14 <- harmonize_cyto(get_col(ucmm_sub, "t(4:14) - FGFR3/IGH fusion at dx"))

t_14_16 <- harmonize_cyto(get_col(ucmm_sub, "hierarchy_corrected_t(14:16) - IGH/MAF translocation at dx"))
if (all(is.na(t_14_16))) t_14_16 <- harmonize_cyto(get_col(ucmm_sub, "t(14:16) - IGH/MAF translocation at dx"))

t_14_20 <- harmonize_cyto(get_col(ucmm_sub, "hierarchy_corrected_t(14:20) - IGH/MAFB translocation at dx"))
if (all(is.na(t_14_20))) t_14_20 <- harmonize_cyto(get_col(ucmm_sub, "t(14:20) - IGH/MAFB translocation at dx"))

del_17p <- harmonize_cyto(get_col(ucmm_sub, "del17p - TP53 Deletion at dx"))
gain_1q <- harmonize_cyto(get_col(ucmm_sub, "gain1q - CKS1B Gain at dx"))
del_1p  <- harmonize_cyto(get_col(ucmm_sub, "del1p - CDKN2C Loss at dx"))

# High-risk cytogenetics: any of t(4;14), t(14;16), t(14;20), del(17p), gain(1q)
# Use HRCA if available
hrca <- to_numeric_safe(get_col(ucmm_sub, "HRCA"))
high_risk_cyto <- ifelse(!is.na(hrca), as.integer(hrca >= 1),
                  ifelse(rowSums(cbind(t_4_14, t_14_16, t_14_20, del_17p, gain_1q),
                                 na.rm = TRUE) > 0, 1, 0))

# --- ASCT ---
has_asct_flag <- as.integer(!is.na(to_days(ucmm_sub$asct_time_abstracted)))

# --- Treatment regimen ---
regimen <- as.character(get_col(ucmm_sub, "regimen"))

# --- Best response ---
response <- as.character(get_col(ucmm_sub, "response_abstracted"))

# --- BMI ---
bmi <- to_numeric_safe(get_col(ucmm_sub, "bmi_dx_emr"))

# ============================================================
# 5. BUILD CLEANED DATA FRAME (primary cohort)
# ============================================================
ucmm_cleaned <- data.frame(
  # ID (synthetic, sequential -- NOT the institutional ERRC ID)
  subject_id     = sprintf("UCMM-%03d", seq_len(nrow(ucmm_sub))),
  source         = "UCMM",

  # Demographics
  age            = round(age_final, 1),
  sex            = sex_clean,
  race           = race_raw,
  ethnicity      = as.character(get_col(ucmm_sub, "eth_composite")),

  # Disease staging
  iss            = iss_final,
  riss           = riss_derived,
  ecog           = NA_character_,   # not available in UCMM

  # Cytogenetics / FISH
  high_risk_cyto = high_risk_cyto,
  t_4_14         = t_4_14,
  t_14_16        = t_14_16,
  t_14_20        = t_14_20,
  del_17p        = del_17p,
  gain_1q        = gain_1q,
  del_1p         = del_1p,
  double_hit     = NA_character_,   # UCMM uses HRCA012; map 2+ as double-hit
  # Override: derive double_hit from HRCA
  # (will set below)

  # Heavy/light chains
  heavy_chain    = NA_character_,   # not directly available in same format
  light_chain    = NA_character_,

  # Treatment details
  transplant_off_protocol = has_asct_flag,
  best_response  = response,

  # Survival outcomes (from tx start date)
  os_months      = round(os_months, 2),
  os_status      = os_event,
  pfs_months     = round(pfs_months, 2),
  pfs_status     = pfs_event,

  # NOTE: calendar dates (c1d1/enrollment/dx/progression/death/last_followup)
  # are intentionally NOT retained here -- only de-identified relative times
  # (os_months, pfs_months) are kept for public release.

  stringsAsFactors = FALSE
)

# Derive double_hit from HRCA (2+ = double-hit)
hrca012 <- to_numeric_safe(get_col(ucmm_sub, "HRCA012"))
ucmm_cleaned$double_hit <- ifelse(!is.na(hrca012) & hrca012 >= 2, "Yes",
                            ifelse(!is.na(hrca012), "No", NA))

# ============================================================
# 6. DEFINE THE TWO ANALYSIS COHORTS
# ============================================================
# PRIMARY   = every patient passing restriction steps 1-4 (any regimen; E-Rd excluded)
# SECONDARY = the primary cohort restricted to KRd or Rd only  (= Cohort 3)
# (Cohort 1 = "Has Tx" n=303 and Cohort 2 = "excl. ASCT on/before" n=254 are the
#  funnel checkpoints shown in the codebook; the analysis uses primary + secondary.)
ucmm_primary   <- ucmm_cleaned
sec_idx        <- which(regimen %in% c("KRd", "Rd"))
ucmm_secondary <- ucmm_cleaned[sec_idx, ]
cat(sprintf("\nPRIMARY   cohort: n=%d (all regimens; E-Rd excluded)\n", nrow(ucmm_primary)))
cat(sprintf("SECONDARY cohort: n=%d (Cohort 3: KRd %d + Rd %d)\n",
            nrow(ucmm_secondary), sum(regimen == "KRd"), sum(regimen == "Rd")))
#    PRIMARY ............................................. n = 253
#    SECONDARY (Cohort 3: KRd 35 + Rd 13) ................ n =  48

# ============================================================
# 7. VALIDATE + SAVE  (filename carries the sample size)
# ============================================================
out_dir <- file.path(mainDir, "data_cleaned")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

save_cohort <- function(df, label, reg) {
  cat(sprintf("\n=== %s cohort: %d subjects x %d columns ===\n",
              label, nrow(df), ncol(df)))
  cat(sprintf("  OS events:  %d / %d (%.0f%%)\n",
              sum(df$os_status == 1), nrow(df), 100 * mean(df$os_status == 1)))
  cat(sprintf("  PFS events: %d / %d (%.0f%%)\n",
              sum(df$pfs_status == 1), nrow(df), 100 * mean(df$pfs_status == 1)))
  cat(sprintf("  Median OS:  %.1f months | Median PFS: %.1f months\n",
              median(df$os_months, na.rm = TRUE), median(df$pfs_months, na.rm = TRUE)))
  cat("  ASCT (transplant):\n"); print(table(df$transplant_off_protocol, useNA = "ifany"))
  cat("  Regimen (top 10):\n"); print(head(sort(table(reg), decreasing = TRUE), 10))
  ucmm_cleaned <- df  # object name kept as ucmm_cleaned for downstream loaders
  out_file <- file.path(out_dir, sprintf("ucmm_cleaned_n%d.RData", nrow(df)))
  save(ucmm_cleaned, file = out_file)
  cat(sprintf("  Saved: %s\n", out_file))
}

save_cohort(ucmm_primary,   "PRIMARY",   regimen)
save_cohort(ucmm_secondary, "SECONDARY", regimen[sec_idx])
