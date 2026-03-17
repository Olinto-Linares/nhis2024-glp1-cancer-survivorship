# ==============================================================================
# NHIS 2024 Analysis Pipeline — Block 00: Data Loading & Sample Construction
# ==============================================================================
# Study: Disparities in GLP-1 Receptor Agonist Use Among Adults With Diabetes
#        and Cancer History: Findings From the 2024 National Health Interview Survey
#
# Description:
#   Loads the NHIS 2024 Sample Adult file and applies the pre-specified
#   exclusion sequence to construct the final analytic sample.
#
# Exclusion sequence:
#   E1. Missing / refused / don't know diabetes status (DIBEV_A not in {1,2})
#   E2. Missing cancer history among those with diagnosed diabetes
#        (CANEV_A not in {1,2})
#   E3. Non-melanoma skin cancer only (SKNNMCAN_A == 1 AND NUMCAN_A == 1)
#   E4. Missing GLP-1 RA outcome (DIBGLP_A not in {1,2})
#        Reported separately by cancer history group to assess
#        differential non-response.
#
# Input:  data/adult24.csv  (NHIS 2024 Sample Adult file)
# Output: final_df  — final analytic dataset passed to Block 01
#
# Data source:
#   Download the NHIS 2024 Adult File from:
#   https://www.cdc.gov/nchs/nhis/2024nhis.htm
#   Place adult24.csv in the data/ folder before running.
#
# Run order: Block_00 → Block_01 → Block_02 → ... → Block_10
# ==============================================================================

# ------------------------------------------------------------------------------
# 0-A  Packages & data path
# ------------------------------------------------------------------------------
library(tidyverse)

if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
library(here)

data_path <- here("data", "adult24.csv")

if (!file.exists(data_path)) {
  stop(
    "NHIS 2024 data file not found.\n",
    "Please download the NHIS 2024 Sample Adult File from:\n",
    "https://www.cdc.gov/nchs/nhis/2024nhis.htm\n",
    "and place adult24.csv in the data/ folder at the project root.\n",
    "See data/README.md for full instructions."
  )
}

df        <- read.csv(data_path, stringsAsFactors = FALSE)
n_total   <- nrow(df)

cat("Full NHIS 2024 adult sample:", n_total, "\n\n")

# ------------------------------------------------------------------------------
# 0-B  E1 — Restrict to adults with diagnosed diabetes (DIBEV_A == 1)
# ------------------------------------------------------------------------------
diabetes_complete     <- df %>% filter(DIBEV_A %in% c(1, 2))
n_e1_missing_diabetes <- n_total - nrow(diabetes_complete)

diabetes_df    <- diabetes_complete %>% filter(DIBEV_A == 1)
n_not_diabetic <- nrow(diabetes_complete) - nrow(diabetes_df)

cat("--- Step 1: Diabetes status ---\n")
cat("  Adults with complete diabetes response:           ", nrow(diabetes_complete), "\n")
cat("  E1 — Excluded (missing diabetes status):          ", n_e1_missing_diabetes, "\n")
cat("  Not diabetic:                                     ", n_not_diabetic, "\n")
cat("  Diagnosed diabetes:                               ", nrow(diabetes_df), "\n\n")

# ------------------------------------------------------------------------------
# 0-C  E2 — Exclude missing cancer history among diabetes patients
# ------------------------------------------------------------------------------
n_diabetes_start        <- nrow(diabetes_df)
diabetes_cancer_complete <- diabetes_df %>% filter(CANEV_A %in% c(1, 2))
n_e2_missing_cancer     <- n_diabetes_start - nrow(diabetes_cancer_complete)

cat("--- Step 2: Cancer history completeness ---\n")
cat("  Diabetes patients entering this step:             ", n_diabetes_start, "\n")
cat("  E2 — Excluded (missing cancer history):           ", n_e2_missing_cancer, "\n")
cat("  Remaining after E2:                               ", nrow(diabetes_cancer_complete), "\n\n")

diabetes_no_cancer  <- diabetes_cancer_complete %>% filter(CANEV_A == 2)
diabetes_any_cancer <- diabetes_cancer_complete %>% filter(CANEV_A == 1)

cat("  No cancer history:                               ", nrow(diabetes_no_cancer), "\n")
cat("  Any cancer history:                              ", nrow(diabetes_any_cancer), "\n\n")

# ------------------------------------------------------------------------------
# 0-D  E3 — Exclude non-melanoma skin cancer only
#           Consistent with NCI SEER cancer survivorship conventions
# ------------------------------------------------------------------------------
nonmel_only <- diabetes_any_cancer %>%
  filter(SKNNMCAN_A == 1 & NUMCAN_A == 1)

n_e3_nonmel <- nrow(nonmel_only)

diabetes_cancer_included <- diabetes_any_cancer %>%
  filter(!(SKNNMCAN_A == 1 & NUMCAN_A == 1))

cat("--- Step 3: Non-melanoma skin cancer only exclusion ---\n")
cat("  Cancer patients entering this step:               ", nrow(diabetes_any_cancer), "\n")
cat("  E3 — Excluded (only non-melanoma skin cancer):   ", n_e3_nonmel, "\n")
cat("  Cancer survivors retained:                        ", nrow(diabetes_cancer_included), "\n\n")

cat("  Retained survivors by number of cancer diagnoses:\n")
print(table(diabetes_cancer_included$NUMCAN_A, useNA = "ifany"))
cat("\n")

# ------------------------------------------------------------------------------
# 0-E  E4 — Combine cohort; exclude missing GLP-1 RA outcome
#           Reported separately by cancer history to assess differential
#           non-response before applying the exclusion.
# ------------------------------------------------------------------------------
cohort_pre_e4          <- bind_rows(diabetes_no_cancer, diabetes_cancer_included)
n_cohort_pre_e4        <- nrow(cohort_pre_e4)
n_cohort_no_cancer_pre <- nrow(diabetes_no_cancer)
n_cohort_cancer_pre    <- nrow(diabetes_cancer_included)

n_e4_no_cancer <- cohort_pre_e4 %>%
  filter(CANEV_A == 2, !DIBGLP_A %in% c(1, 2)) %>% nrow()

n_e4_cancer <- cohort_pre_e4 %>%
  filter(CANEV_A == 1, !DIBGLP_A %in% c(1, 2)) %>% nrow()

n_e4_total <- n_e4_no_cancer + n_e4_cancer

cat("--- Step 4: GLP-1 RA outcome completeness (E4) ---\n")
cat("  Combined cohort entering this step:               ", n_cohort_pre_e4, "\n")
cat("    No cancer history:                              ", n_cohort_no_cancer_pre, "\n")
cat("    Cancer history:                                 ", n_cohort_cancer_pre, "\n\n")
cat("  E4 — Excluded (missing GLP-1 RA outcome):\n")
cat("    No cancer history:                              ", n_e4_no_cancer,
    sprintf("(%.1f%% of no-cancer group)\n", 100 * n_e4_no_cancer / n_cohort_no_cancer_pre))
cat("    Cancer history:                                 ", n_e4_cancer,
    sprintf("(%.1f%% of cancer group)\n", 100 * n_e4_cancer / n_cohort_cancer_pre))
cat("    Total excluded (E4):                            ", n_e4_total, "\n\n")

final_df <- cohort_pre_e4 %>% filter(DIBGLP_A %in% c(1, 2))

n_final_no_cancer <- sum(final_df$CANEV_A == 2)
n_final_cancer    <- sum(final_df$CANEV_A == 1)

cat("  Final analytic sample after E4:                  ", nrow(final_df), "\n")
cat("    No cancer history:                              ", n_final_no_cancer,
    sprintf("(%.1f%%)\n", 100 * n_final_no_cancer / nrow(final_df)))
cat("    Cancer history:                                 ", n_final_cancer,
    sprintf("(%.1f%%)\n", 100 * n_final_cancer / nrow(final_df)))
cat("\n")

# ------------------------------------------------------------------------------
# 0-F  Exclusion flow summary
# ------------------------------------------------------------------------------
cat("=============================================================================\n")
cat("EXCLUSION FLOW SUMMARY\n")
cat("=============================================================================\n")
cat(sprintf("  Full NHIS 2024 adult sample:                      %d\n", n_total))
cat(sprintf("  E1  Missing / invalid diabetes status:           -%d\n", n_e1_missing_diabetes))
cat(sprintf("      Not diabetic (out of scope):                 -%d\n", n_not_diabetic))
cat(sprintf("      -> Diagnosed diabetes:                        %d\n", n_diabetes_start))
cat(sprintf("  E2  Missing cancer history:                      -%d\n", n_e2_missing_cancer))
cat(sprintf("      -> Diabetes + complete cancer info:           %d\n", nrow(diabetes_cancer_complete)))
cat(sprintf("  E3  Non-melanoma skin cancer only:               -%d\n", n_e3_nonmel))
cat(sprintf("      -> Combined cohort (pre-outcome):             %d\n", n_cohort_pre_e4))
cat(sprintf("          No cancer history:                        %d\n", n_cohort_no_cancer_pre))
cat(sprintf("          Cancer history:                           %d\n", n_cohort_cancer_pre))
cat(sprintf("  E4  Missing GLP-1 RA outcome:                    -%d\n", n_e4_total))
cat(sprintf("          No cancer history:                        %d\n", n_e4_no_cancer))
cat(sprintf("          Cancer history:                           %d\n", n_e4_cancer))
cat(sprintf("  -------------------------------------------------------\n"))
cat(sprintf("  FINAL ANALYTIC SAMPLE:                            %d\n", nrow(final_df)))
cat(sprintf("    No cancer history:                              %d (%.1f%%)\n",
            n_final_no_cancer, 100 * n_final_no_cancer / nrow(final_df)))
cat(sprintf("    Cancer history (survivors):                     %d (%.1f%%)\n",
            n_final_cancer, 100 * n_final_cancer / nrow(final_df)))
cat("=============================================================================\n\n")

# ------------------------------------------------------------------------------
# 0-G  Verification checks
# ------------------------------------------------------------------------------
cat("--- Verification checks (all should be 0) ---\n")

cat("  Non-melanoma-only cases in final_df:             ",
    sum(final_df$SKNNMCAN_A == 1 & final_df$NUMCAN_A == 1, na.rm = TRUE), "\n")
cat("  Missing GLP-1 outcome in final_df:               ",
    sum(!final_df$DIBGLP_A %in% c(1, 2)), "\n")
cat("  Missing DIBEV_A in final_df:                     ",
    sum(!final_df$DIBEV_A %in% c(1, 2), na.rm = TRUE), "\n")
cat("  Missing CANEV_A in final_df:                     ",
    sum(!final_df$CANEV_A %in% c(1, 2), na.rm = TRUE), "\n")

cat("\n=============================================================================\n")
cat("Block 00 complete — final_df ready (n =", nrow(final_df), ")\n")
cat("=============================================================================\n")
