rm(list = ls())

mainDir <- "/Users/oliviazhang/Desktop/mapbart-historical-borrowing/mapbart-case-study-mm"
clean_dir <- file.path(mainDir, "data_cleaned")

# ============================================================
# 1. LOAD CLEANED DATASETS
# ============================================================
elokrd_files <- list.files(clean_dir, pattern = "^elokrd_cleaned_n.*\\.RData$",
                           full.names = TRUE)
if (length(elokrd_files) == 0) stop("EloKRD cleaned file not found. Run data_cleaning_elokrd.R first.")
elokrd_file <- elokrd_files[which.max(file.mtime(elokrd_files))]
load(elokrd_file)  # loads 'elokrd'
cat(sprintf("EloKRD: %d x %d\n", nrow(elokrd), ncol(elokrd)))

# ============================================================
# Merge EloKRD against EACH cleaned UCMM file found in data_cleaned/.
# Cohorts are identified purely by file name (ucmm_cleaned_n<N>.RData);
# each produces merged_elokrd_ucmm_n<N>.RData.
# ============================================================
ucmm_files <- list.files(clean_dir, pattern = "^ucmm_cleaned_n[0-9]+\\.RData$",
                         full.names = TRUE)
if (length(ucmm_files) == 0)
  stop("No ucmm_cleaned_n*.RData in ", clean_dir, "  (run data_cleaning_ucmm.R first)")

for (ucmm_file in ucmm_files) {
cat(sprintf("\n################ MERGE: %s ################\n", basename(ucmm_file)))

.e <- new.env(); load(ucmm_file, envir = .e)   # loads 'ucmm_cleaned'
ucmm_cleaned <- .e$ucmm_cleaned
cat(sprintf("UCMM (%s): %d x %d\n", basename(ucmm_file),
            nrow(ucmm_cleaned), ncol(ucmm_cleaned)))

# ============================================================
# 2. HARMONIZE FACTOR LEVELS BEFORE FILTERING
# ============================================================

# --- Sex ---
harmonize_sex <- function(x) {
  x <- trimws(as.character(x))
  ifelse(grepl("^[Ff]", x), "Female",
  ifelse(grepl("^[Mm]", x), "Male", NA))
}
elokrd$sex       <- harmonize_sex(elokrd$sex)
ucmm_cleaned$sex <- harmonize_sex(ucmm_cleaned$sex)

# --- Race ---
harmonize_race <- function(x) {
  x <- trimws(tolower(as.character(x)))
  ifelse(grepl("white|caucasian", x), "White",
  ifelse(grepl("black|african", x), "Black",
  ifelse(grepl("asian", x), "Asian",
  ifelse(grepl("hispanic|latino", x), "Hispanic",
  ifelse(x %in% c("", "na", "unknown", "not reported", "other"), "Other/Unknown",
         "Other/Unknown")))))
}
elokrd$race       <- harmonize_race(elokrd$race)
ucmm_cleaned$race <- harmonize_race(ucmm_cleaned$race)

# --- ISS ---
harmonize_iss <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^0-9IViv]", "", x)
  ifelse(grepl("^1$|^I$", x), "I",
  ifelse(grepl("^2$|^II$", x), "II",
  ifelse(grepl("^3$|^III$", x), "III", NA)))
}
elokrd$iss       <- harmonize_iss(elokrd$iss)
ucmm_cleaned$iss <- harmonize_iss(ucmm_cleaned$iss)

# --- R-ISS ---
harmonize_riss <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[^0-9IViv]", "", x)
  ifelse(grepl("^1$|^I$", x), "I",
  ifelse(grepl("^2$|^II$", x), "II",
  ifelse(grepl("^3$|^III$", x), "III", NA)))
}
elokrd$riss       <- harmonize_riss(elokrd$riss)
ucmm_cleaned$riss <- harmonize_riss(ucmm_cleaned$riss)

# --- ECOG ---
elokrd$ecog       <- as.character(elokrd$ecog)
ucmm_cleaned$ecog <- as.character(ucmm_cleaned$ecog)

# --- Double-hit ---
elokrd$double_hit <- ifelse(grepl("[Yy]es|^1$", as.character(elokrd$double_hit)), "Yes",
                     ifelse(grepl("[Nn]o|^0$", as.character(elokrd$double_hit)), "No", NA))

# --- Best response ---
harmonize_response <- function(x) {
  x <- trimws(toupper(as.character(x)))
  ifelse(grepl("^SCR|STRINGENT", x), "sCR",
  ifelse(grepl("^CR$|COMPLETE RESP", x), "CR",
  ifelse(grepl("^VGPR|VERY GOOD", x), "VGPR",
  ifelse(grepl("^PR$|PARTIAL RESP", x), "PR",
  ifelse(grepl("^SD|STABLE", x), "SD",
  ifelse(grepl("^PD|PROGRESS", x), "PD", NA))))))
}
elokrd$best_response       <- harmonize_response(elokrd$best_response)
ucmm_cleaned$best_response <- harmonize_response(ucmm_cleaned$best_response)

# ============================================================
# 3. IDENTIFY VARIABLES WITH 100% COMPLETENESS IN EloKRD
# ============================================================
common_cols <- intersect(names(elokrd), names(ucmm_cleaned))

# Exclude ID/metadata and date columns (always kept or handled separately)
id_meta_cols <- c("subject_id", "source")
date_cols <- c("c1d1_date", "enrollment_date", "dx_date",
               "progression_date", "death_date", "last_followup")

candidate_cols <- setdiff(common_cols, c(id_meta_cols, date_cols))

n_elokrd <- nrow(elokrd)
cat(sprintf("\n=== STEP 1: Variable completeness in EloKRD (n=%d) ===\n", n_elokrd))
elokrd_complete <- c()
elokrd_drop <- c()
for (v in candidate_cols) {
  n_complete <- sum(!is.na(elokrd[[v]]))
  status <- ifelse(n_complete == n_elokrd, "KEEP", "DROP")
  cat(sprintf("  %-30s %2d/%2d  %s\n", v, n_complete, n_elokrd, status))
  if (n_complete == n_elokrd) {
    elokrd_complete <- c(elokrd_complete, v)
  } else {
    elokrd_drop <- c(elokrd_drop, v)
  }
}
cat(sprintf("\n%d/%d in EloKRD: %d variables\n", n_elokrd, n_elokrd, length(elokrd_complete)))
cat(sprintf("Dropped (incomplete in EloKRD): %s\n", paste(elokrd_drop, collapse = ", ")))

# ============================================================
# 4. FROM THE 30/30 VARS, FIND THE LARGEST SUBSET THAT
#    MAXIMIZES UCMM COMPLETE CASES
#    Strategy: iteratively drop the variable with the most
#    missingness in UCMM until we reach a good tradeoff,
#    then require complete cases on what remains.
# ============================================================
cat("\n=== STEP 2: UCMM completeness on 100% EloKRD variables ===\n")
for (v in elokrd_complete) {
  n_avail <- sum(!is.na(ucmm_cleaned[[v]]))
  cat(sprintf("  %-30s %3d/%3d (%.0f%%)\n", v, n_avail, nrow(ucmm_cleaned),
              100 * n_avail / nrow(ucmm_cleaned)))
}

# Iteratively drop worst UCMM-coverage variable and show how many
# UCMM complete cases we'd get at each step
cat("\n=== STEP 3: Iterative variable trimming (maximize UCMM n) ===\n")
cat(sprintf("  %-4s  %-30s  %8s  %7s  %s\n",
            "Step", "Dropped variable", "UCMM avail", "UCMM n", "Variables remaining"))

remaining <- elokrd_complete
step <- 0

# Show baseline
n_cc <- sum(complete.cases(ucmm_cleaned[, remaining]))
cat(sprintf("  %4d  %-30s  %8s  %7d  %d vars\n",
            step, "(none - all 30/30 vars)", "", n_cc, length(remaining)))

# Iteratively drop the variable with the lowest UCMM availability
# until complete-case n stabilizes
results <- list()
results[[1]] <- list(dropped = "(none)", vars = remaining, n_cc = n_cc)

while (length(remaining) > 0 && n_cc < nrow(ucmm_cleaned)) {
  # Find variable with fewest non-NA in UCMM
  avail <- sapply(remaining, function(v) sum(!is.na(ucmm_cleaned[[v]])))
  worst <- names(which.min(avail))
  worst_avail <- avail[worst]

  # Drop it
  remaining <- setdiff(remaining, worst)
  step <- step + 1

  if (length(remaining) == 0) break

  n_cc <- sum(complete.cases(ucmm_cleaned[, remaining]))
  cat(sprintf("  %4d  %-30s  %3d/%3d    %7d  %d vars\n",
              step, worst, worst_avail, nrow(ucmm_cleaned), n_cc, length(remaining)))

  results[[step + 1]] <- list(dropped = worst, vars = remaining, n_cc = n_cc)
}

# Pick the step that gives the best tradeoff: maximize
# (number of UCMM patients) while keeping as many variables as possible
# Rule: pick the first step where UCMM n >= 30 (matching EloKRD size minimum)
# and we keep the most covariates.
# If multiple steps have the same n, prefer more variables.
cat("\n=== STEP 4: Select optimal variable set ===\n")

# Find steps where we have a reasonable UCMM n
# We want to retain as many UCMM patients as possible while keeping covariates
min_ucmm_n <- nrow(ucmm_cleaned)  # use all available UCMM patients as target
best_idx <- NULL
for (i in seq_along(results)) {
  r <- results[[i]]
  # Find the first result where all UCMM patients are complete cases
  if (r$n_cc >= min_ucmm_n) {
    best_idx <- i
    break
  }
}

# If nothing >= 30, take the one with max n_cc
if (is.null(best_idx)) {
  ncc_vals <- sapply(results, function(r) r$n_cc)
  best_idx <- which.max(ncc_vals)
}

best <- results[[best_idx]]
keep_vars <- best$vars

n_cc_selected <- best$n_cc
cat(sprintf("Selected: %d variables, %d UCMM complete cases\n",
            length(keep_vars), n_cc_selected))
cat(sprintf("Variables kept: %s\n", paste(keep_vars, collapse = ", ")))

# Now subset UCMM
ucmm_complete_flag <- complete.cases(ucmm_cleaned[, keep_vars])
n_before <- nrow(ucmm_cleaned)
ucmm_cleaned <- ucmm_cleaned[ucmm_complete_flag, ]
n_after <- nrow(ucmm_cleaned)

cat(sprintf("\nUCMM: %d -> %d (dropped %d with missingness)\n",
            n_before, n_after, n_before - n_after))

# Verify: no missingness remains
cat("\n=== Verification: completeness after subsetting ===\n")
for (v in keep_vars) {
  n_avail <- sum(!is.na(ucmm_cleaned[[v]]))
  cat(sprintf("  %-30s %3d/%3d\n", v, n_avail, nrow(ucmm_cleaned)))
}

# ============================================================
# 5. MERGE (row-bind) — kept variables + ID columns only
# ============================================================
merge_cols <- c(id_meta_cols, keep_vars)

merged <- rbind(elokrd[, merge_cols], ucmm_cleaned[, merge_cols])

# Add treatment group indicator (1 = EloKRD trial, 0 = UCMM control)
merged$trt <- ifelse(merged$source == "EloKRD", 1, 0)

cat(sprintf("\n=== MERGED DATASET: %d subjects (%d EloKRD + %d UCMM) x %d columns ===\n",
            nrow(merged), sum(merged$trt == 1), sum(merged$trt == 0), ncol(merged)))

# ============================================================
# 6. SUMMARY STATISTICS
# ============================================================
cat("\n--- Demographics by group ---\n")
for (grp in c("EloKRD", "UCMM")) {
  d <- merged[merged$source == grp, ]
  cat(sprintf("\n  %s (n=%d):\n", grp, nrow(d)))
  cat(sprintf("    Age: median %.1f (range %.0f-%.0f)\n",
              median(d$age, na.rm = TRUE), min(d$age, na.rm = TRUE), max(d$age, na.rm = TRUE)))
  cat(sprintf("    Male: %d (%.0f%%)\n",
              sum(d$sex == "Male", na.rm = TRUE),
              100 * mean(d$sex == "Male", na.rm = TRUE)))
  if ("iss" %in% names(d)) {
    cat("    ISS: "); print(table(d$iss, useNA = "ifany"))
  }
  if ("high_risk_cyto" %in% names(d)) {
    cat(sprintf("    High-risk cyto: %d (%.0f%%)\n",
                sum(d$high_risk_cyto == 1, na.rm = TRUE),
                100 * mean(d$high_risk_cyto == 1, na.rm = TRUE)))
  }
  cat(sprintf("    OS events: %d (%.0f%%)\n",
              sum(d$os_status == 1, na.rm = TRUE),
              100 * mean(d$os_status == 1, na.rm = TRUE)))
  cat(sprintf("    Median OS: %.1f months\n", median(d$os_months, na.rm = TRUE)))
  cat(sprintf("    PFS events: %d (%.0f%%)\n",
              sum(d$pfs_status == 1, na.rm = TRUE),
              100 * mean(d$pfs_status == 1, na.rm = TRUE)))
  cat(sprintf("    Median PFS: %.1f months\n", median(d$pfs_months, na.rm = TRUE)))
}

# Final completeness check — should be 100% everywhere
cat("\n--- Final variable completeness (should be 100%) ---\n")
for (v in keep_vars) {
  n_avail_e <- sum(!is.na(merged[[v]][merged$trt == 1]))
  n_avail_u <- sum(!is.na(merged[[v]][merged$trt == 0]))
  cat(sprintf("  %-30s  EloKRD: %2d/%2d  UCMM: %3d/%3d\n",
              v, n_avail_e, sum(merged$trt == 1), n_avail_u, sum(merged$trt == 0)))
}

# ============================================================
# 7. SAVE MERGED DATASET
# ============================================================
# Filename carries the total sample size, which identifies the cohort.
out_file <- file.path(clean_dir, sprintf("merged_elokrd_ucmm_n%d.RData", nrow(merged)))
save(merged, file = out_file)
cat(sprintf("\nSaved: %s\n", out_file))

cat("\n=== Column listing ===\n")
cat(paste(sprintf("  [%2d] %s", seq_along(names(merged)), names(merged)), collapse = "\n"), "\n")

cat("\n=== DONE ===\n")
cat("Merged dataset ready for synthetic control analysis.\n")
cat("Key columns:\n")
cat("  trt         : 1 = EloKRD (trial), 0 = UCMM (control)\n")
cat("  os_months   : OS time from tx start (C1D1 for EloKRD, induction for UCMM)\n")
cat("  os_status   : 1 = death, 0 = censored\n")
cat("  pfs_months  : PFS time from tx start\n")
cat("  pfs_status  : 1 = progression or death, 0 = censored\n")

}  # end loop over UCMM cleaned files
