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



# ==============================================================================
# NHIS 2024 Analysis Pipeline — Block 01: Variable Construction
# ==============================================================================
# Input:  final_df  (from Block 00)
# Output: analytic_df  — final_df with all analysis variables appended
#
# Variables constructed:
#   Outcome       : glp1_use
#   Exposure      : cancer_history
#   Demographics  : age_years, age_group, sex, race_ethnicity, region
#   SES           : education, income_fpl_2cat (primary), income_fpl_4cat (S1)
#   Insurance     : insurance_type  (3 categories)
#   Clinical      : bmi_category, ascvd, ckd, hypertension, hyperlipidemia
#   Access        : outpatient_12m  (descriptive only; excluded from models)
#   Exploratory   : cost_barrier, spd_k6  (Model C only)
#   Diabetes type : DIBTYPE_A  (retained for sensitivity analysis S6)
#
# Pre-specified analytic adaptations (documented for transparency):
#   A1. Insurance: planned 4-category -> 3 categories (Private/Medicare/Medicaid).
#       Uninsured set to NA and excluded from regression models.
#       Rationale: only 6 uninsured cancer survivors; insufficient for stable
#       variance estimation.
#   A2. BMI: planned 4-category (<25, 25-29.9, 30-34.9, >=35) ->
#       3 categories (<30, 30-34.9, >=35).
#       Rationale: sparse cell at BMI <25 within the cancer survivor stratum.
#   A3. Race/ethnicity: planned 5-category -> 4 categories (NH White, NH Black,
#       Hispanic, Other). Non-Hispanic Asian collapsed into "Other" due to
#       sparse cells in the cancer survivor stratum. Retained in descriptive
#       analyses; collapsed for regression.
#   A4. outpatient_12m excluded from regression models (98.5% prevalence;
#       insufficient discriminatory variation). Included in Table 1.
# ==============================================================================

library(dplyr)

if (!exists("final_df")) stop("final_df not found. Run Block 00 first.")

analytic_df <- final_df

cat("=============================================================================\n")
cat("Block 01: Variable Construction\n")
cat("Input sample size:", nrow(analytic_df), "\n")
cat("=============================================================================\n\n")

# ==============================================================================
# 1-1  EXPOSURE — Cancer history
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    cancer_history = factor(
      case_when(
        CANEV_A == 1 ~ "Cancer history",
        CANEV_A == 2 ~ "No cancer history",
        TRUE         ~ NA_character_
      ),
      levels = c("No cancer history", "Cancer history")
    )
  )

cat("1-1  Cancer history (exposure)\n")
print(table(analytic_df$cancer_history, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-2  OUTCOME — GLP-1 RA use (DIBGLP_A)
#      NHIS 2024 variable DIBGLP_A — new Emerging Content in the 2024 cycle.
#      Question: "[Other than insulin,] are you NOW taking any injectable
#      medications to lower your blood sugar or lose weight?"
#      1 = Yes, 2 = No
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    glp1_use = factor(
      case_when(
        DIBGLP_A == 1 ~ "Yes",
        DIBGLP_A == 2 ~ "No",
        TRUE          ~ NA_character_
      ),
      levels = c("No", "Yes")
    )
  )

cat("1-2  GLP-1 RA use (outcome)\n")
print(table(analytic_df$glp1_use, useNA = "ifany"))
cat("\nBy cancer history:\n")
print(addmargins(table(analytic_df$glp1_use, analytic_df$cancer_history, useNA = "ifany")))
cat("\n")

# ==============================================================================
# 1-3  DEMOGRAPHICS — Age (AGEP_A)
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    age_years = AGEP_A,
    age_group = factor(
      case_when(
        AGEP_A >= 18 & AGEP_A <= 44 ~ "18\u201344",
        AGEP_A >= 45 & AGEP_A <= 64 ~ "45\u201364",
        AGEP_A >= 65               ~ "65+",
        TRUE                       ~ NA_character_
      ),
      levels = c("18\u201344", "45\u201364", "65+")
    )
  )

cat("1-3  Age\n")
cat("  Mean (SD):", round(mean(analytic_df$age_years, na.rm = TRUE), 1),
    "(", round(sd(analytic_df$age_years, na.rm = TRUE), 1), ")\n")
print(table(analytic_df$age_group, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-4  DEMOGRAPHICS — Sex (SEX_A)
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    sex = factor(
      case_when(
        SEX_A == 1 ~ "Male",
        SEX_A == 2 ~ "Female",
        TRUE       ~ NA_character_
      ),
      levels = c("Male", "Female")
    )
  )

cat("1-4  Sex\n")
print(table(analytic_df$sex, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-5  DEMOGRAPHICS — Race/ethnicity (HISPALLP_A)
#      Analytic adaptation A3: non-Hispanic Asian and other groups collapsed
#      into "Other" due to sparse cells in the cancer survivor stratum.
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    race_ethnicity = factor(
      case_when(
        HISPALLP_A == 1               ~ "Hispanic",
        HISPALLP_A == 2               ~ "Non-Hispanic White",
        HISPALLP_A == 3               ~ "Non-Hispanic Black",
        HISPALLP_A %in% c(4,5,6,7,8) ~ "Other",
        TRUE                           ~ NA_character_
      ),
      levels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other")
    )
  )

cat("1-5  Race/ethnicity (4 categories; see analytic adaptation A3)\n")
print(table(analytic_df$race_ethnicity, useNA = "ifany"))
cat("\nBy cancer history:\n")
print(addmargins(table(analytic_df$race_ethnicity, analytic_df$cancer_history, useNA = "ifany")))
cat("\n")

# ==============================================================================
# 1-6  DEMOGRAPHICS — Geographic region (REGION)
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    region = factor(
      case_when(
        REGION == 1 ~ "Northeast",
        REGION == 2 ~ "Midwest",
        REGION == 3 ~ "South",
        REGION == 4 ~ "West",
        TRUE        ~ NA_character_
      ),
      levels = c("Northeast", "Midwest", "South", "West")
    )
  )

cat("1-6  Geographic region\n")
print(table(analytic_df$region, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-7  SES — Education (EDUCP_A)
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    education = factor(
      case_when(
        EDUCP_A %in% c(1, 2, 3) ~ "High school or less",
        EDUCP_A %in% c(4, 5)    ~ "Some college",
        EDUCP_A %in% c(6, 7, 8) ~ "Bachelor's+",
        TRUE                    ~ NA_character_
      ),
      levels = c("High school or less", "Some college", "Bachelor's+")
    )
  )

cat("1-7  Education\n")
print(table(analytic_df$education, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-8  SES — Household income as % FPL (POVRATTC_A)
#      Primary specification (Models A-B): binary (<200% vs >=200%)
#      Sensitivity S1: 4-category breakdown
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    income_fpl_2cat = factor(
      case_when(
        is.na(POVRATTC_A) | POVRATTC_A < 0 ~ NA_character_,
        POVRATTC_A < 2.0                    ~ "<200% FPL",
        POVRATTC_A >= 2.0                   ~ "\u2265200% FPL",
        TRUE                                ~ NA_character_
      ),
      levels = c("<200% FPL", "\u2265200% FPL")
    ),
    income_fpl_4cat = factor(
      case_when(
        is.na(POVRATTC_A) | POVRATTC_A < 0 ~ NA_character_,
        POVRATTC_A < 1.0                    ~ "<100% FPL",
        POVRATTC_A < 2.0                    ~ "100\u2013199% FPL",
        POVRATTC_A < 4.0                    ~ "200\u2013399% FPL",
        POVRATTC_A >= 4.0                   ~ "\u2265400% FPL",
        TRUE                                ~ NA_character_
      ),
      levels = c("<100% FPL", "100\u2013199% FPL", "200\u2013399% FPL", "\u2265400% FPL")
    )
  )

cat("1-8  Income (% FPL)\n")
cat("  Primary (2-category):\n")
print(table(analytic_df$income_fpl_2cat, useNA = "ifany"))
cat("  Sensitivity S1 (4-category):\n")
print(table(analytic_df$income_fpl_4cat, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-9  INSURANCE — Insurance type (PRIVATE_A, MEDICARE_A, MEDICAID_A, NOTCOV_A)
#      Analytic adaptation A1: uninsured set to NA and excluded from regression.
#      n = 182 uninsured overall; n = 6 in cancer survivor stratum.
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    insurance_type = factor(
      case_when(
        PRIVATE_A  == 1 ~ "Private",
        MEDICARE_A == 1 ~ "Medicare",
        MEDICAID_A == 1 ~ "Medicaid",
        NOTCOV_A   == 1 ~ NA_character_,
        TRUE            ~ NA_character_
      ),
      levels = c("Private", "Medicare", "Medicaid")
    )
  )

n_uninsured <- sum(analytic_df$NOTCOV_A == 1, na.rm = TRUE)

cat("1-9  Insurance type (3-category; see analytic adaptation A1)\n")
print(table(analytic_df$insurance_type, useNA = "ifany"))
cat("\n  Uninsured (set to NA, excluded from regression):", n_uninsured, "\n\n")

# ==============================================================================
# 1-10  CLINICAL — BMI category (BMICAT_A)
#       Analytic adaptation A2: <25 and 25-29.9 collapsed to <30
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    bmi_category = factor(
      case_when(
        BMICAT_A %in% c(1, 2) ~ "<30",
        BMICAT_A == 3         ~ "30\u201334.9",
        BMICAT_A == 4         ~ "\u226535",
        TRUE                  ~ NA_character_
      ),
      levels = c("<30", "30\u201334.9", "\u226535")
    )
  )

cat("1-10  BMI category (3-category; see analytic adaptation A2)\n")
print(table(analytic_df$bmi_category, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-11  CLINICAL — ASCVD composite
#       Components: CHD (CHDEV_A), angina (ANGEV_A), MI (MIEV_A), stroke (STREV_A)
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    ascvd = factor(
      case_when(
        CHDEV_A == 1 | ANGEV_A == 1 | MIEV_A == 1 | STREV_A == 1 ~ "Yes",
        CHDEV_A == 2 & ANGEV_A == 2 & MIEV_A == 2 & STREV_A == 2 ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    )
  )

cat("1-11  ASCVD composite (CHD / angina / MI / stroke)\n")
print(table(analytic_df$ascvd, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-12  CLINICAL — Chronic kidney disease (KIDWEAKEV_A)
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    ckd = factor(
      case_when(
        KIDWEAKEV_A == 1 ~ "Yes",
        KIDWEAKEV_A == 2 ~ "No",
        TRUE             ~ NA_character_
      ),
      levels = c("No", "Yes")
    )
  )

cat("1-12  Chronic kidney disease\n")
print(table(analytic_df$ckd, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-13  CLINICAL — Hypertension (HYPEV_A)  [Table 1 / Sensitivity S2]
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    hypertension = factor(
      case_when(
        HYPEV_A == 1 ~ "Yes",
        HYPEV_A == 2 ~ "No",
        TRUE         ~ NA_character_
      ),
      levels = c("No", "Yes")
    )
  )

cat("1-13  Hypertension\n")
print(table(analytic_df$hypertension, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-14  CLINICAL — Hyperlipidemia (CHLEV_A)  [Table 1 / Sensitivity S2]
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    hyperlipidemia = factor(
      case_when(
        CHLEV_A == 1 ~ "Yes",
        CHLEV_A == 2 ~ "No",
        TRUE         ~ NA_character_
      ),
      levels = c("No", "Yes")
    )
  )

cat("1-14  Hyperlipidemia\n")
print(table(analytic_df$hyperlipidemia, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-15  ACCESS — Outpatient visit in past 12 months (LASTDR_A)
#       Analytic adaptation A4: excluded from regression (98.5% prevalence);
#       included in Table 1 as descriptive access indicator.
#       1-2 = visited within past year ("Yes"); 3-5 = >1 year or never ("No")
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    outpatient_12m = factor(
      case_when(
        LASTDR_A %in% c(1, 2)   ~ "Yes",
        LASTDR_A %in% c(3, 4, 5) ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    )
  )

cat("1-15  Outpatient visit past 12 months (descriptive only; see adaptation A4)\n")
print(table(analytic_df$outpatient_12m, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-16  EXPLORATORY (Model C) — Cost barrier and serious psychological distress
#       Cost barrier: any delayed (MEDDL12M_A) or unmet (MEDNG12M_A) care due
#       to cost in the past 12 months.
#       Serious psychological distress: Kessler-6 scale (K6SPD_A).
# ==============================================================================
analytic_df <- analytic_df %>%
  mutate(
    cost_barrier = factor(
      case_when(
        MEDDL12M_A == 1 | MEDNG12M_A == 1 ~ "Yes",
        MEDDL12M_A == 2 & MEDNG12M_A == 2 ~ "No",
        TRUE ~ NA_character_
      ),
      levels = c("No", "Yes")
    ),
    spd_k6 = factor(
      case_when(
        K6SPD_A == 1 ~ "Yes",
        K6SPD_A == 2 ~ "No",
        TRUE         ~ NA_character_
      ),
      levels = c("No", "Yes")
    )
  )

cat("1-16  Exploratory variables (Model C)\n")
cat("  Cost barrier:\n")
print(table(analytic_df$cost_barrier, useNA = "ifany"))
cat("  Serious psychological distress (Kessler-6):\n")
print(table(analytic_df$spd_k6, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-17  Retain diabetes type for sensitivity analysis S6
# ==============================================================================
analytic_df$DIBTYPE_A <- final_df$DIBTYPE_A

cat("1-17  Diabetes type (retained for sensitivity analysis S6)\n")
print(table(analytic_df$DIBTYPE_A, useNA = "ifany"))
cat("\n")

# ==============================================================================
# 1-18  Missing data summary
# ==============================================================================
analytic_vars <- c(
  "glp1_use", "cancer_history",
  "age_group", "sex", "race_ethnicity", "region",
  "education", "income_fpl_2cat", "income_fpl_4cat",
  "insurance_type",
  "bmi_category", "ascvd", "ckd", "hypertension", "hyperlipidemia",
  "outpatient_12m", "cost_barrier", "spd_k6"
)

missing_summary <- sapply(analytic_vars, function(v) sum(is.na(analytic_df[[v]])))

cat("=============================================================================\n")
cat("Missing data summary (n missing per variable)\n")
cat("=============================================================================\n")
print(missing_summary)
cat("\n")

primary_vars <- c(
  "glp1_use", "cancer_history", "age_years", "sex", "race_ethnicity",
  "education", "income_fpl_2cat", "insurance_type",
  "bmi_category", "ascvd", "ckd"
)

analytic_df$complete_primary <- complete.cases(analytic_df[, primary_vars])
n_complete_primary <- sum(analytic_df$complete_primary)

cat("Complete cases on primary analytic variables (Models A/B):", n_complete_primary,
    sprintf("(%.1f%% of %d)\n", 100 * n_complete_primary / nrow(analytic_df), nrow(analytic_df)))

cat("\n=============================================================================\n")
cat("Block 01 complete — analytic_df ready (n =", nrow(analytic_df), ")\n")
cat("=============================================================================\n\n")
cat("VARIABLE INVENTORY\n")
cat("  Outcome        : glp1_use\n")
cat("  Exposure       : cancer_history\n")
cat("  Demographics   : age_years, age_group, sex, race_ethnicity, region\n")
cat("  SES            : education, income_fpl_2cat [primary], income_fpl_4cat [S1]\n")
cat("  Insurance      : insurance_type [3-category; see adaptation A1]\n")
cat("  Clinical       : bmi_category [3-category; A2], ascvd, ckd\n")
cat("  Table 1 extras : hypertension, hyperlipidemia, outpatient_12m [A4]\n")
cat("  Exploratory    : cost_barrier, spd_k6 [Model C]\n")
cat("  Diabetes type  : DIBTYPE_A [retained for sensitivity S6]\n\n")
cat("ANALYTIC ADAPTATIONS\n")
cat("  A1. Insurance 4 -> 3 categories (uninsured n=6 in cancer stratum)\n")
cat("  A2. BMI 4 -> 3 categories (<25 + 25-29.9 collapsed to <30)\n")
cat("  A3. Race/ethnicity 5 -> 4 categories (NH Asian collapsed into Other)\n")
cat("  A4. outpatient_12m excluded from regression (98.5% prevalence)\n")
cat("=============================================================================\n")

# ==============================================================================
# NHIS 2024 Analysis Pipeline — Block 02: Survey Design Setup
# ==============================================================================
# Input:  analytic_df (from Block 01)
# Output: nhis_design  — svydesign object used in all subsequent blocks
#
# The NHIS uses a multistage probability sampling design with stratification
# and clustering. All analyses must account for this design using:
#   - PSUs (PPSU): primary sampling units
#   - Strata (PSTRAT): design strata
#   - Annual weights (WTFA_A): full-year survey weights
#
# Uninsured respondents (NA in insurance_type) remain in nhis_design but are
# automatically excluded from models via complete-case analysis.
# ==============================================================================

library(survey)

if (!exists("analytic_df")) stop("analytic_df not found. Run Block 01 first.")

cat("=============================================================================\n")
cat("Block 02: Survey Design Setup\n")
cat("=============================================================================\n\n")

nhis_design <- svydesign(
  ids     = ~PPSU,
  strata  = ~PSTRAT,
  weights = ~WTFA_A,
  nest    = TRUE,
  data    = analytic_df
)

cat("Survey design created\n\n")
cat("Unweighted n:               ", nrow(analytic_df), "\n")
cat("Weighted population estimate:", format(round(sum(weights(nhis_design))), big.mark = ","), "\n\n")

cat("Cancer history — unweighted counts:\n")
print(table(analytic_df$cancer_history, useNA = "ifany"))
cat("\nCancer history — weighted estimates:\n")
print(svytable(~cancer_history, nhis_design))
cat("\nUninsured (NA in insurance_type, excluded from regression):",
    sum(is.na(analytic_df$insurance_type)), "\n")

cat("\n=============================================================================\n")
cat("Block 02 complete — nhis_design ready\n")
cat("=============================================================================\n")


# ==============================================================================
# NHIS 2024 Analysis Pipeline — Block 03: Table 1 — Baseline Characteristics
# ==============================================================================
# Input:  analytic_df (Block 01), nhis_design (Block 02)
# Output: table1_baseline_characteristics.docx
#
# Notes:
#   - income_fpl_4cat used for descriptives (4-category breakdown)
#   - income_fpl_2cat used in regression models (Block 04)
#   - outpatient_12m included as descriptive-only access indicator
#   - Uninsured excluded from insurance rows (NA in insurance_type)
# ==============================================================================

library(survey)
library(dplyr)
library(flextable)
library(officer)

if (!exists("nhis_design")) stop("nhis_design not found. Run Block 02 first.")

cat("=============================================================================\n")
cat("Block 03: Table 1 — Baseline Characteristics\n")
cat("=============================================================================\n\n")

n_total     <- nrow(analytic_df)
n_cancer    <- sum(analytic_df$cancer_history == "Cancer history",    na.rm = TRUE)
n_no_cancer <- sum(analytic_df$cancer_history == "No cancer history", na.rm = TRUE)

cat("Unweighted sample sizes:\n")
cat("  Total:            ", n_total, "\n")
cat("  Cancer history:   ", n_cancer, "\n")
cat("  No cancer history:", n_no_cancer, "\n\n")

# ==============================================================================
# Helper functions
# ==============================================================================
fmt_pct_se  <- function(est, se) { if (any(is.na(c(est,se)))) return("—"); sprintf("%.1f (%.1f)", est*100, se*100) }
fmt_mean_se <- function(est, se) { if (any(is.na(c(est,se)))) return("—"); sprintf("%.1f (%.1f)", est, se) }
fmt_p       <- function(p) { if (is.na(p)) return("—"); if (p < 0.001) return("<.001"); sprintf("%.3f", p) }

pull_est <- function(by_res, group, varname, level) {
  col <- paste0(varname, level)
  row <- by_res[by_res$cancer_history == group, , drop = FALSE]
  val <- suppressWarnings(as.numeric(row[[col]]))
  if (!is.na(val)) return(val)
  pat  <- paste0("^", gsub("([.^$*+?{}\\[\\]|()\\/ ])", "\\\\\1", col), "$")
  cand <- grep(pat, names(row), value = TRUE)
  cand <- cand[!startsWith(cand, "se.")]
  if (length(cand) == 1) return(as.numeric(row[[cand]]))
  NA_real_
}
pull_se <- function(by_res, group, varname, level) {
  col  <- paste0("se.", varname, level)
  row  <- by_res[by_res$cancer_history == group, , drop = FALSE]
  val  <- suppressWarnings(as.numeric(row[[col]]))
  if (!is.na(val)) return(val)
  pat  <- paste0("^", gsub("([.^$*+?{}\\[\\]|()\\/ ])", "\\\\\1", col), "$")
  cand <- grep(pat, names(row), value = TRUE)
  if (length(cand) == 1) return(as.numeric(row[[cand]]))
  NA_real_
}
pull_ov_est <- function(res, varname, level) {
  col <- paste0(varname, level)
  val <- suppressWarnings(as.numeric(coef(res)[col]))
  if (!is.na(val)) return(val)
  pat  <- paste0("^", gsub("([.^$*+?{}\\[\\]|()\\/ ])", "\\\\\1", col), "$")
  cand <- grep(pat, names(coef(res)), value = TRUE)
  if (length(cand) == 1) return(as.numeric(coef(res)[cand]))
  NA_real_
}
pull_ov_se <- function(res, varname, level) {
  col <- paste0(varname, level)
  val <- suppressWarnings(as.numeric(SE(res)[col]))
  if (!is.na(val)) return(val)
  pat  <- paste0("^", gsub("([.^$*+?{}\\[\\]|()\\/ ])", "\\\\\1", col), "$")
  cand <- grep(pat, names(SE(res)), value = TRUE)
  if (length(cand) == 1) return(as.numeric(SE(res)[cand]))
  NA_real_
}

make_cat_rows <- function(char_label, varname, levels_vec, display_vec, by_res, ov_res, p_val) {
  bind_rows(lapply(seq_along(levels_vec), function(i) {
    lv <- levels_vec[i]
    data.frame(
      Characteristic = if (i == 1) char_label else "",
      Category       = display_vec[i],
      Overall        = fmt_pct_se(pull_ov_est(ov_res, varname, lv), pull_ov_se(ov_res, varname, lv)),
      Cancer_history = fmt_pct_se(pull_est(by_res, "Cancer history",    varname, lv),
                                  pull_se(by_res,  "Cancer history",    varname, lv)),
      No_cancer      = fmt_pct_se(pull_est(by_res, "No cancer history", varname, lv),
                                  pull_se(by_res,  "No cancer history", varname, lv)),
      P_value        = if (i == 1) p_val else "",
      stringsAsFactors = FALSE
    )
  }))
}

make_bin_row <- function(char_label, varname, by_res, ov_res, p_val) {
  data.frame(
    Characteristic = char_label, Category = "Yes",
    Overall        = fmt_pct_se(pull_ov_est(ov_res, varname, "Yes"), pull_ov_se(ov_res, varname, "Yes")),
    Cancer_history = fmt_pct_se(pull_est(by_res, "Cancer history",    varname, "Yes"),
                                pull_se(by_res,  "Cancer history",    varname, "Yes")),
    No_cancer      = fmt_pct_se(pull_est(by_res, "No cancer history", varname, "Yes"),
                                pull_se(by_res,  "No cancer history", varname, "Yes")),
    P_value = p_val, stringsAsFactors = FALSE
  )
}

make_header_row <- function(label) {
  data.frame(Characteristic = label, Category = "", Overall = "",
             Cancer_history = "", No_cancer = "", P_value = "", stringsAsFactors = FALSE)
}

# ==============================================================================
# Compute weighted statistics
# ==============================================================================
cat("Computing survey-weighted statistics...\n\n")

by_grp <- function(var) svyby(as.formula(paste0("~", var)), ~cancer_history, nhis_design, svymean, na.rm = TRUE, vartype = "se")
ov_mean <- function(var) svymean(as.formula(paste0("~", var)), nhis_design, na.rm = TRUE)
chi_p   <- function(var) fmt_p(svychisq(as.formula(paste0("~", var, "+cancer_history")), nhis_design, statistic = "Chisq")$p.value)
t_p     <- function(var) fmt_p(svyttest(as.formula(paste0(var, "~cancer_history")), nhis_design)$p.value)

age_by      <- svyby(~age_years, ~cancer_history, nhis_design, svymean, na.rm = TRUE)
age_ov      <- svymean(~age_years, nhis_design, na.rm = TRUE)
age_grp_by  <- by_grp("age_group");      age_grp_ov  <- ov_mean("age_group")
sex_by      <- by_grp("sex");            sex_ov      <- ov_mean("sex")
race_by     <- by_grp("race_ethnicity"); race_ov     <- ov_mean("race_ethnicity")
edu_by      <- by_grp("education");      edu_ov      <- ov_mean("education")
inc_by      <- by_grp("income_fpl_4cat");inc_ov      <- ov_mean("income_fpl_4cat")
ins_by      <- by_grp("insurance_type"); ins_ov      <- ov_mean("insurance_type")
reg_by      <- by_grp("region");         reg_ov      <- ov_mean("region")
bmi_by      <- by_grp("bmi_category");   bmi_ov      <- ov_mean("bmi_category")
asc_by      <- by_grp("ascvd");          asc_ov      <- ov_mean("ascvd")
ckd_by      <- by_grp("ckd");            ckd_ov      <- ov_mean("ckd")
htn_by      <- by_grp("hypertension");   htn_ov      <- ov_mean("hypertension")
hld_by      <- by_grp("hyperlipidemia"); hld_ov      <- ov_mean("hyperlipidemia")
op_by       <- by_grp("outpatient_12m"); op_ov       <- ov_mean("outpatient_12m")
cost_by     <- by_grp("cost_barrier");   cost_ov     <- ov_mean("cost_barrier")
spd_by      <- by_grp("spd_k6");         spd_ov      <- ov_mean("spd_k6")
glp1_by     <- by_grp("glp1_use");       glp1_ov     <- ov_mean("glp1_use")

cat("Statistics computed\n\n")

# ==============================================================================
# Build Table 1 row by row
# ==============================================================================
rows <- list()

rows[[1]] <- data.frame(
  Characteristic = "Sample size (unweighted n)", Category = "—",
  Overall = as.character(n_total), Cancer_history = as.character(n_cancer),
  No_cancer = as.character(n_no_cancer), P_value = "—", stringsAsFactors = FALSE
)

rows[[2]]  <- make_header_row("Demographics")
ca_age     <- age_by[age_by$cancer_history == "Cancer history", ]
nc_age     <- age_by[age_by$cancer_history == "No cancer history", ]
rows[[3]]  <- data.frame(
  Characteristic = "Age, years", Category = "Mean (SE)",
  Overall        = fmt_mean_se(coef(age_ov), SE(age_ov)),
  Cancer_history = fmt_mean_se(ca_age$age_years, ca_age$se),
  No_cancer      = fmt_mean_se(nc_age$age_years, nc_age$se),
  P_value        = t_p("age_years"), stringsAsFactors = FALSE
)
rows[[4]]  <- make_cat_rows("Age group", "age_group",
                            c("18\u201344", "45\u201364", "65+"), c("18\u201344", "45\u201364", "\u226565"),
                            age_grp_by, age_grp_ov, chi_p("age_group"))
rows[[5]]  <- make_cat_rows("Sex", "sex", c("Male", "Female"), c("Male", "Female"),
                            sex_by, sex_ov, chi_p("sex"))
rows[[6]]  <- make_cat_rows("Race/ethnicity", "race_ethnicity",
                            c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other"),
                            c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other\u2020"),
                            race_by, race_ov, chi_p("race_ethnicity"))

rows[[7]]  <- make_header_row("Socioeconomic factors")
rows[[8]]  <- make_cat_rows("Education", "education",
                            c("High school or less", "Some college", "Bachelor's+"),
                            c("High school or less", "Some college", "Bachelor's+"),
                            edu_by, edu_ov, chi_p("education"))
rows[[9]]  <- make_cat_rows("Income (% FPL)", "income_fpl_4cat",
                            c("<100% FPL", "100\u2013199% FPL", "200\u2013399% FPL", "\u2265400% FPL"),
                            c("<100%", "100\u2013199%", "200\u2013399%", "\u2265400%"),
                            inc_by, inc_ov, chi_p("income_fpl_4cat"))

rows[[10]] <- make_header_row("Insurance")
rows[[11]] <- make_cat_rows("Insurance type\u2021", "insurance_type",
                            c("Private", "Medicare", "Medicaid"), c("Private", "Medicare", "Medicaid"),
                            ins_by, ins_ov, chi_p("insurance_type"))
rows[[12]] <- make_cat_rows("Region", "region",
                            c("Northeast", "Midwest", "South", "West"), c("Northeast", "Midwest", "South", "West"),
                            reg_by, reg_ov, chi_p("region"))

rows[[13]] <- make_header_row("Cardiometabolic comorbidities")
rows[[14]] <- make_cat_rows("BMI category (kg/m\u00b2)", "bmi_category",
                            c("<30", "30\u201334.9", "\u226535"), c("<30", "30\u201334.9", "\u226535"),
                            bmi_by, bmi_ov, chi_p("bmi_category"))
rows[[15]] <- make_bin_row("ASCVD (CHD/MI or stroke)", "ascvd",  asc_by, asc_ov, chi_p("ascvd"))
rows[[16]] <- make_bin_row("Chronic kidney disease",    "ckd",    ckd_by, ckd_ov, chi_p("ckd"))
rows[[17]] <- make_bin_row("Hypertension",              "hypertension", htn_by, htn_ov, chi_p("hypertension"))
rows[[18]] <- make_bin_row("Hyperlipidemia",            "hyperlipidemia", hld_by, hld_ov, chi_p("hyperlipidemia"))

rows[[19]] <- make_header_row("Access indicators (descriptive)")
rows[[20]] <- make_bin_row("Outpatient visit past 12 months", "outpatient_12m",
                           op_by, op_ov, chi_p("outpatient_12m"))
rows[[21]] <- make_bin_row("Cost barrier (unmet/delayed care)", "cost_barrier",
                           cost_by, cost_ov, chi_p("cost_barrier"))
rows[[22]] <- make_bin_row("Serious psychological distress (K6)", "spd_k6",
                           spd_by, spd_ov, chi_p("spd_k6"))

rows[[23]] <- make_header_row("Outcome")
rows[[24]] <- make_bin_row("Current GLP-1 RA use", "glp1_use",
                           glp1_by, glp1_ov, chi_p("glp1_use"))

table1_df <- bind_rows(rows)
colnames(table1_df) <- c("Characteristic", "Category", "Overall",
                         "Cancer history", "No cancer history", "P value")
cat("Table 1 built:", nrow(table1_df), "rows\n\n")

# ==============================================================================
# Export to Word
# ==============================================================================
header_rows <- which(table1_df$Characteristic %in% c(
  "Demographics", "Socioeconomic factors", "Insurance",
  "Cardiometabolic comorbidities", "Access indicators (descriptive)", "Outcome"
))

ft1 <- flextable(table1_df) %>%
  theme_booktabs() %>% autofit() %>%
  align(align = "left",   j = 1:2, part = "all") %>%
  align(align = "center", j = 3:6, part = "body") %>%
  align(align = "center", j = 3:6, part = "header") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  bold(part = "header") %>%
  bold(i = header_rows, j = 1) %>%
  bg(i = header_rows, bg = "#F2F2F2") %>%
  set_header_labels(
    Characteristic    = "Characteristic",
    Category          = "Category",
    Overall           = paste0("Overall\n(n = ", n_total, ")"),
    `Cancer history`  = paste0("Cancer history\n(n = ", n_cancer, ")"),
    `No cancer history` = paste0("No cancer history\n(n = ", n_no_cancer, ")"),
    `P value`         = "P value"
  )

footnote_t1 <- paste(
  "Cells are weighted % (SE) unless otherwise noted; continuous variables shown as weighted mean (SE).",
  "P values from survey-weighted \u03c7\u00b2 tests (categorical) and survey-weighted t-tests (continuous).",
  "Cancer history: any cancer diagnosis excluding nonmelanoma skin cancer only.",
  "\u2020 \'Other\' race/ethnicity includes non-Hispanic Asian, non-Hispanic Pacific Islander,",
  "American Indian or Alaska Native, other single race, and multiracial respondents;",
  "collapsed due to sparse cells in the cancer survivor stratum.",
  "\u2021 Uninsured respondents (n =", sum(is.na(analytic_df$insurance_type)),
  ") excluded from insurance rows and regression models.",
  "BMI category collapsed to 3 groups due to sparse cells.",
  "Cost barrier: any delayed or unmet care due to cost in past 12 months.",
  "Serious psychological distress: Kessler-6 scale.",
  "SE = standard error; FPL = federal poverty level; ASCVD = atherosclerotic cardiovascular disease;",
  "GLP-1 RA = glucagon-like peptide-1 receptor agonist."
)

doc1 <- read_docx() %>%
  body_add_par(
    "Table 1. Weighted Baseline Characteristics of U.S. Adults With Diagnosed Diabetes, by Cancer History (NHIS 2024)",
    style = "heading 1"
  ) %>%
  body_add_par("") %>%
  body_add_flextable(ft1) %>%
  body_add_par("") %>%
  body_add_par(footnote_t1, style = "Normal")

print(doc1, target = "table1_baseline_characteristics.docx")

cat("\n=============================================================================\n")
cat("Block 03 complete — table1_baseline_characteristics.docx saved\n")
cat("=============================================================================\n")

# ==============================================================================
# NHIS 2024 Analysis Pipeline — Block 04: Table 2 — Primary Analysis
# ==============================================================================
# Input:  analytic_df (Block 01), nhis_design (Block 02)
# Output: table2_glp1_models.docx
#         model_a, model_b  (used in Blocks 05-10)
#
# Models:
#   Unadjusted : GLP-1 RA ~ cancer_history only
#   Model A    : + age + sex + race/ethnicity + education + income +
#                  insurance + region  (composition-standardized)
#   Model B    : Model A + BMI + ASCVD + CKD  (need-standardized)
#
# Marginal standardization (G-computation) produces adjusted predicted
# probabilities averaged over the complete-case covariate distribution,
# fully survey-weighted via Taylor series linearization.
# Absolute difference 95% CI uses the delta method (conservative).
# ==============================================================================

library(survey)
library(dplyr)
library(flextable)
library(officer)

if (!exists("nhis_design")) stop("nhis_design not found. Run Block 02 first.")

cat("=============================================================================\n")
cat("Block 04: Table 2 — GLP-1 RA Use Analysis\n")
cat("=============================================================================\n\n")

fmt_est_ci <- function(est, lo, hi) {
  if (any(is.na(c(est, lo, hi)))) return("\u2014")
  sprintf("%.1f (%.1f, %.1f)", est, lo, hi)
}
fmt_or_ci <- function(or, lo, hi) {
  if (any(is.na(c(or, lo, hi)))) return("\u2014")
  sprintf("%.2f (%.2f, %.2f)", or, lo, hi)
}
fmt_pp_ci <- function(diff, lo, hi) {
  if (any(is.na(c(diff, lo, hi)))) return("\u2014")
  sprintf("%.1f (%.1f, %.1f)", diff, lo, hi)
}

# Sample sizes
n_cancer    <- sum(analytic_df$cancer_history == "Cancer history",    na.rm = TRUE)
n_no_cancer <- sum(analytic_df$cancer_history == "No cancer history", na.rm = TRUE)
n_total     <- nrow(analytic_df)
wt_dist     <- prop.table(svytable(~cancer_history, nhis_design)) * 100
pct_cancer  <- as.numeric(wt_dist["Cancer history"])
pct_no_cancer <- as.numeric(wt_dist["No cancer history"])

cat("  Total:", n_total, "| Cancer:", n_cancer,
    sprintf("(%.1f%% weighted)", pct_cancer),
    "| No cancer:", n_no_cancer,
    sprintf("(%.1f%% weighted)\n\n", pct_no_cancer))

# Unadjusted prevalence
glp1_by_grp <- svyby(~glp1_use, ~cancer_history, nhis_design, svymean, na.rm = TRUE)
get_prev <- function(group) {
  row <- glp1_by_grp[glp1_by_grp$cancer_history == group, ]
  est <- 100 * as.numeric(row[["glp1_useYes"]])
  se  <- 100 * as.numeric(row[["se.glp1_useYes"]])
  list(est = est, lo = est - 1.96*se, hi = est + 1.96*se, se = se)
}
unadj_cancer    <- get_prev("Cancer history")
unadj_no_cancer <- get_prev("No cancer history")
cat("Unadjusted prevalence:\n")
cat("  Cancer history:   ", fmt_est_ci(unadj_cancer$est,    unadj_cancer$lo,    unadj_cancer$hi), "\n")
cat("  No cancer history:", fmt_est_ci(unadj_no_cancer$est, unadj_no_cancer$lo, unadj_no_cancer$hi), "\n\n")

# Complete-case flags injected into design$variables
vars_a <- c("glp1_use","cancer_history","age_years","sex","race_ethnicity",
            "education","income_fpl_2cat","insurance_type","region")
vars_b <- c(vars_a, "bmi_category","ascvd","ckd")
cc_a_flag <- complete.cases(analytic_df[, vars_a])
cc_b_flag <- complete.cases(analytic_df[, vars_b])
nhis_design$variables$cc_a <- cc_a_flag
nhis_design$variables$cc_b <- cc_b_flag
cat("Complete cases: Model A n =", sum(cc_a_flag), "| Model B n =", sum(cc_b_flag), "\n\n")

# Model A
model_a <- svyglm(
  glp1_use ~ cancer_history + age_years + sex + race_ethnicity +
    education + income_fpl_2cat + insurance_type + region,
  design = nhis_design, family = quasibinomial()
)

design_cc_a <- subset(nhis_design, cc_a)
cc_a_df     <- design_cc_a$variables
design_cc_a$variables$pred_cancer_a <- as.numeric(predict(model_a,
                                                          newdata = cc_a_df %>% mutate(cancer_history = factor("Cancer history",    levels = levels(analytic_df$cancer_history))),
                                                          type = "response"))
design_cc_a$variables$pred_no_cancer_a <- as.numeric(predict(model_a,
                                                             newdata = cc_a_df %>% mutate(cancer_history = factor("No cancer history", levels = levels(analytic_df$cancer_history))),
                                                             type = "response"))
std_a_ca <- list(est = 100*as.numeric(coef(svymean(~pred_cancer_a,    design_cc_a, na.rm=TRUE))),
                 se  = 100*as.numeric(SE(svymean(~pred_cancer_a,      design_cc_a, na.rm=TRUE))))
std_a_nc <- list(est = 100*as.numeric(coef(svymean(~pred_no_cancer_a, design_cc_a, na.rm=TRUE))),
                 se  = 100*as.numeric(SE(svymean(~pred_no_cancer_a,   design_cc_a, na.rm=TRUE))))
std_a_ca <- within(std_a_ca, { lo <- est - 1.96*se; hi <- est + 1.96*se })
std_a_nc <- within(std_a_nc, { lo <- est - 1.96*se; hi <- est + 1.96*se })
cat("Model A standardized: Cancer", fmt_est_ci(std_a_ca$est, std_a_ca$lo, std_a_ca$hi),
    "| No cancer", fmt_est_ci(std_a_nc$est, std_a_nc$lo, std_a_nc$hi), "\n")

# Model B
model_b <- svyglm(
  glp1_use ~ cancer_history + age_years + sex + race_ethnicity +
    education + income_fpl_2cat + insurance_type + region +
    bmi_category + ascvd + ckd,
  design = nhis_design, family = quasibinomial()
)

design_cc_b <- subset(nhis_design, cc_b)
cc_b_df     <- design_cc_b$variables
design_cc_b$variables$pred_cancer_b <- as.numeric(predict(model_b,
                                                          newdata = cc_b_df %>% mutate(cancer_history = factor("Cancer history",    levels = levels(analytic_df$cancer_history))),
                                                          type = "response"))
design_cc_b$variables$pred_no_cancer_b <- as.numeric(predict(model_b,
                                                             newdata = cc_b_df %>% mutate(cancer_history = factor("No cancer history", levels = levels(analytic_df$cancer_history))),
                                                             type = "response"))
std_b_ca <- list(est = 100*as.numeric(coef(svymean(~pred_cancer_b,    design_cc_b, na.rm=TRUE))),
                 se  = 100*as.numeric(SE(svymean(~pred_cancer_b,      design_cc_b, na.rm=TRUE))))
std_b_nc <- list(est = 100*as.numeric(coef(svymean(~pred_no_cancer_b, design_cc_b, na.rm=TRUE))),
                 se  = 100*as.numeric(SE(svymean(~pred_no_cancer_b,   design_cc_b, na.rm=TRUE))))
std_b_ca <- within(std_b_ca, { lo <- est - 1.96*se; hi <- est + 1.96*se })
std_b_nc <- within(std_b_nc, { lo <- est - 1.96*se; hi <- est + 1.96*se })
cat("Model B standardized: Cancer", fmt_est_ci(std_b_ca$est, std_b_ca$lo, std_b_ca$hi),
    "| No cancer", fmt_est_ci(std_b_nc$est, std_b_nc$lo, std_b_nc$hi), "\n")

# Absolute difference (delta method)
diff_b    <- std_b_ca$est - std_b_nc$est
se_diff_b <- sqrt(std_b_ca$se^2 + std_b_nc$se^2)
diff_lo   <- diff_b - 1.96 * se_diff_b
diff_hi   <- diff_b + 1.96 * se_diff_b
cat("Absolute difference:", fmt_pp_ci(diff_b, diff_lo, diff_hi), "pp\n")

# Adjusted OR
coef_ca <- coef(model_b)["cancer_historyCancer history"]
se_ca   <- sqrt(vcov(model_b)["cancer_historyCancer history","cancer_historyCancer history"])
or_b    <- exp(coef_ca); or_lo <- exp(coef_ca - 1.96*se_ca); or_hi <- exp(coef_ca + 1.96*se_ca)
p_b     <- summary(model_b)$coefficients["cancer_historyCancer history","Pr(>|t|)"]
cat("Adjusted OR:", fmt_or_ci(or_b, or_lo, or_hi), sprintf("P = %.3f\n\n", p_b))

# Build Table 2
table2_df <- data.frame(
  Group    = c("Cancer history", "No cancer history"),
  n_wpct   = c(sprintf("%d (%.1f%%)", n_cancer, pct_cancer),
               sprintf("%d (%.1f%%)", n_no_cancer, pct_no_cancer)),
  Unadj    = c(fmt_est_ci(unadj_cancer$est,    unadj_cancer$lo,    unadj_cancer$hi),
               fmt_est_ci(unadj_no_cancer$est, unadj_no_cancer$lo, unadj_no_cancer$hi)),
  Model_A  = c(fmt_est_ci(std_a_ca$est, std_a_ca$lo, std_a_ca$hi),
               fmt_est_ci(std_a_nc$est, std_a_nc$lo, std_a_nc$hi)),
  Model_B  = c(fmt_est_ci(std_b_ca$est, std_b_ca$lo, std_b_ca$hi),
               fmt_est_ci(std_b_nc$est, std_b_nc$lo, std_b_nc$hi)),
  Abs_diff = c("\u2014", fmt_pp_ci(diff_b, diff_lo, diff_hi)),
  OR_B     = c("\u2014", fmt_or_ci(or_b, or_lo, or_hi)),
  stringsAsFactors = FALSE
)
colnames(table2_df) <- c(
  "Group", "n (weighted %)",
  "Unadjusted prevalence\n% (95% CI)",
  "Model A\nstandardized % (95% CI)",
  "Model B\nstandardized % (95% CI)",
  "Absolute difference\n(Cancer \u2212 No cancer), pp (95% CI)",
  "Adjusted OR\n(Model B) (95% CI)"
)

ft2 <- flextable(table2_df) %>%
  theme_booktabs() %>% autofit() %>%
  align(align = "left",   j = 1,   part = "all") %>%
  align(align = "center", j = 2:7, part = "body") %>%
  align(align = "center", j = 2:7, part = "header") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  bold(part = "header") %>% bold(j = 1)

footnote_t2 <- paste(
  "Survey-weighted estimates accounting for NHIS 2024 complex sampling design (PSUs, strata, annual weights).",
  "Model A: adjusted for age (continuous), sex, race/ethnicity (4 categories), education (3 categories),",
  "income (<200% vs \u2265200% FPL), insurance type (Private/Medicare/Medicaid), and region.",
  paste0("Uninsured respondents (n = ", sum(is.na(analytic_df$insurance_type)),
         ") excluded from Models A and B via complete-case analysis."),
  "Model B: Model A + BMI category (3 groups: <30, 30\u201334.9, \u226535 kg/m\u00b2), ASCVD (CHD/MI or stroke), and CKD.",
  "Standardized prevalence computed via marginal standardization, survey-weighted.",
  "Absolute difference 95% CI computed via delta method (independence approximation; conservative).",
  "Adjusted OR: Model B logistic regression coefficient for cancer history (Cancer vs No cancer).",
  "Model B comorbidities may be influenced by cancer treatment history (prespecified overadjustment caveat).",
  "CI = confidence interval; pp = percentage points; OR = odds ratio; FPL = federal poverty level."
)

doc2 <- read_docx() %>%
  body_add_par("Table 2. GLP-1 Receptor Agonist Use Among U.S. Adults With Diagnosed Diabetes by Cancer History: Unadjusted and Standardized Estimates (NHIS 2024)", style = "heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft2) %>%
  body_add_par("") %>%
  body_add_par(footnote_t2, style = "Normal")

print(doc2, target = "table2_glp1_models.docx")

cat("\n=============================================================================\n")
cat(sprintf("Block 04 complete | Gap: %.1f pp (%.1f, %.1f) | aOR: %.2f (%.2f, %.2f)\n",
            diff_b, diff_lo, diff_hi, or_b, or_lo, or_hi))
cat("table2_glp1_models.docx saved\n")
cat("=============================================================================\n")


# ==============================================================================
# NHIS 2024 Analysis Pipeline — Block 05: Effect Modification Analysis (Table 3)
# ==============================================================================
# Input:  analytic_df (Block 01), nhis_design (Block 02), model_b (Block 04)
# Output: table3_effect_modification.docx
#
# Tests effect modification of the cancer-GLP-1 RA association by:
#   Panel A: Insurance type (Private, Medicare, Medicaid)
#   Panel B: Race/ethnicity (NH White, NH Black, Hispanic; "Other" excluded
#            from stratification due to heterogeneity and small cell size)
#
# Interaction tests use Wald tests (regTermTest) at alpha = 0.10.
# Stratified estimates use marginal standardization (Model B main effects).
# Cell suppression threshold: n >= 30 per cancer history group per stratum,
# consistent with NCHS analytic guidelines for NHIS public-use data.
# ==============================================================================

library(dplyr)
library(survey)
library(flextable)
library(officer)

if (!exists("nhis_design")) stop("nhis_design not found. Run Block 02 first.")
if (!exists("model_b"))     stop("model_b not found. Run Block 04 first.")

NCHS_THRESHOLD <- 30   # Minimum unweighted n per cancer history group per stratum

cat("=============================================================================\n")
cat("Block 05: Effect Modification Analysis\n")
cat(sprintf("Cell suppression threshold: n >= %d per cancer history group\n", NCHS_THRESHOLD))
cat("Race/ethnicity stratification: NH White, NH Black, Hispanic\n")
cat("=============================================================================\n\n")

# ==============================================================================
# 5-0  Formatting helpers
# ==============================================================================
fmt_prev  <- function(p, se) {
  if (any(is.na(c(p, se))) || any(is.nan(c(p, se)))) return("\u2014")
  sprintf("%.1f (%.1f)", p, se)
}
fmt_diff  <- function(d, lo, hi) {
  if (any(is.na(c(d, lo, hi))) || any(is.nan(c(d, lo, hi)))) return("\u2014")
  sprintf("%.1f (%.1f, %.1f)", d, lo, hi)
}
fmt_p_val <- function(p) {
  if (is.null(p) || is.na(p)) return("\u2014")
  if (p < 0.001) return("<.001")
  sprintf("%.3f", p)
}

# ==============================================================================
# 5-1  Marginal standardization helper
# ==============================================================================
marg_std_subset <- function(model, subset_design, label) {
  model_vars <- c("glp1_use", "cancer_history", "age_years", "sex",
                  "race_ethnicity", "education", "income_fpl_2cat",
                  "insurance_type", "region", "bmi_category", "ascvd", "ckd")
  model_vars_present <- intersect(model_vars, names(subset_design$variables))
  cc_flag <- complete.cases(subset_design$variables[, model_vars_present])
  
  subset_design$variables$cc_flag <- cc_flag
  design_cc <- subset(subset_design, cc_flag)
  cc_df     <- design_cc$variables
  
  n_ca <- sum(cc_df$cancer_history == "Cancer history",    na.rm = TRUE)
  n_nc <- sum(cc_df$cancer_history == "No cancer history", na.rm = TRUE)
  cat(sprintf("    Complete cases: %d (Cancer = %d, No cancer = %d)\n",
              nrow(cc_df), n_ca, n_nc))
  
  if (n_ca < 10 || n_nc < 10) {
    cat("    Too few complete cases — skipping\n")
    return(NULL)
  }
  
  cc_as_ca <- cc_df %>%
    mutate(cancer_history = factor("Cancer history",    levels = levels(analytic_df$cancer_history)))
  cc_as_nc <- cc_df %>%
    mutate(cancer_history = factor("No cancer history", levels = levels(analytic_df$cancer_history)))
  
  design_cc$variables$pred_ca <- as.numeric(predict(model, newdata = cc_as_ca, type = "response"))
  design_cc$variables$pred_nc <- as.numeric(predict(model, newdata = cc_as_nc, type = "response"))
  
  obj_ca <- svymean(~pred_ca, design_cc, na.rm = TRUE)
  obj_nc <- svymean(~pred_nc, design_cc, na.rm = TRUE)
  
  p_ca <- 100 * as.numeric(coef(obj_ca));  se_ca <- 100 * as.numeric(SE(obj_ca))
  p_nc <- 100 * as.numeric(coef(obj_nc));  se_nc <- 100 * as.numeric(SE(obj_nc))
  diff    <- p_ca - p_nc
  se_diff <- sqrt(se_ca^2 + se_nc^2)
  
  list(
    prev_ca = p_ca,  se_ca   = se_ca,
    prev_nc = p_nc,  se_nc   = se_nc,
    diff    = diff,  se_diff = se_diff,
    ci_lo   = diff - 1.96 * se_diff,
    ci_hi   = diff + 1.96 * se_diff,
    n_ca    = n_ca,  n_nc    = n_nc
  )
}

# ==============================================================================
# 5-2  Interaction tests
# ==============================================================================
cat("--- Step 1: Interaction tests ---\n\n")

p_int_ins  <- 1.0
p_int_race <- 1.0

cat("Cancer x Insurance interaction:\n")
tryCatch({
  m_int_ins <- svyglm(
    glp1_use ~ cancer_history * insurance_type +
      age_years + sex + race_ethnicity +
      education + income_fpl_2cat + region +
      bmi_category + ascvd + ckd,
    design = nhis_design, family = quasibinomial()
  )
  tryCatch({
    wt <- regTermTest(m_int_ins, ~cancer_history:insurance_type, method = "Wald")
    p_int_ins <- as.numeric(wt$p)
    cat("  Wald test P for interaction =", fmt_p_val(p_int_ins), "\n\n")
  }, error = function(e) {
    coefs <- summary(m_int_ins)$coefficients
    rows  <- grep("cancer_historyCancer history:insurance_type", rownames(coefs))
    if (length(rows) > 0) p_int_ins <<- min(coefs[rows, "Pr(>|t|)"], na.rm = TRUE)
    cat("  P for interaction (min individual term) =", fmt_p_val(p_int_ins), "\n\n")
  })
}, error = function(e) {
  cat("  Model error:", conditionMessage(e), "\n  P set to 1.0\n\n")
})

cat("Cancer x Race/Ethnicity interaction:\n")
tryCatch({
  d_3race <- subset(nhis_design,
                    race_ethnicity %in% c("Non-Hispanic White",
                                          "Non-Hispanic Black", "Hispanic"))
  m_int_race <- svyglm(
    glp1_use ~ cancer_history * race_ethnicity +
      age_years + sex +
      education + income_fpl_2cat + insurance_type + region +
      bmi_category + ascvd + ckd,
    design = d_3race, family = quasibinomial()
  )
  tryCatch({
    wt <- regTermTest(m_int_race, ~cancer_history:race_ethnicity, method = "Wald")
    p_int_race <- as.numeric(wt$p)
    cat("  Wald test P for interaction =", fmt_p_val(p_int_race), "\n\n")
  }, error = function(e) {
    coefs <- summary(m_int_race)$coefficients
    rows  <- grep("cancer_historyCancer history:race_ethnicity", rownames(coefs))
    if (length(rows) > 0) p_int_race <<- min(coefs[rows, "Pr(>|t|)"], na.rm = TRUE)
    cat("  P for interaction =", fmt_p_val(p_int_race), "\n\n")
  })
}, error = function(e) {
  cat("  Model error:", conditionMessage(e), "\n  P set to 1.0\n\n")
})

p_int_ins  <- ifelse(is.na(p_int_ins),  1.0, p_int_ins)
p_int_race <- ifelse(is.na(p_int_race), 1.0, p_int_race)

THRESHOLD    <- 0.10
include_main <- (p_int_ins < THRESHOLD) | (p_int_race < THRESHOLD)

cat("--- Interaction results ---\n")
cat("  Cancer x Insurance:      P =", fmt_p_val(p_int_ins), "\n")
cat("  Cancer x Race/Ethnicity: P =", fmt_p_val(p_int_race), "(3 groups)\n")
cat("  Table 3 in main text:   ", include_main, "(p < 0.10 criterion)\n\n")

# ==============================================================================
# 5-3  Stratified analysis — by insurance type
# ==============================================================================
cat("--- Step 2: Stratified analysis by insurance type ---\n\n")

insurance_levels <- c("Private", "Medicare", "Medicaid")
ins_results      <- list()

for (ins in insurance_levels) {
  cat("Insurance:", ins, "\n")
  d_ins <- subset(nhis_design, insurance_type == ins)
  n_tot <- nrow(d_ins$variables)
  n_ca  <- sum(d_ins$variables$cancer_history == "Cancer history",    na.rm = TRUE)
  n_nc  <- sum(d_ins$variables$cancer_history == "No cancer history", na.rm = TRUE)
  cat(sprintf("  n total = %d | Cancer = %d | No cancer = %d\n", n_tot, n_ca, n_nc))
  
  if (n_ca < NCHS_THRESHOLD || n_nc < NCHS_THRESHOLD) {
    cat(sprintf("  Suppressed (n < %d per group per NCHS guidelines)\n\n", NCHS_THRESHOLD))
    ins_results[[ins]] <- list(n_tot = n_tot, n_ca = n_ca, n_nc = n_nc, skip = TRUE)
    next
  }
  
  tryCatch({
    m_ins <- svyglm(
      glp1_use ~ cancer_history +
        age_years + sex + race_ethnicity +
        education + income_fpl_2cat + region +
        bmi_category + ascvd + ckd,
      design = d_ins, family = quasibinomial()
    )
    std <- marg_std_subset(m_ins, d_ins, ins)
    if (is.null(std)) {
      ins_results[[ins]] <- list(n_tot = n_tot, n_ca = n_ca, n_nc = n_nc, skip = TRUE)
    } else {
      cat(sprintf("  Cancer: %.1f%% (SE %.1f) | No cancer: %.1f%% (SE %.1f) | Gap: %.1f pp\n",
                  std$prev_ca, std$se_ca, std$prev_nc, std$se_nc, std$diff))
      ins_results[[ins]] <- c(list(n_tot = n_tot, n_ca = n_ca, n_nc = n_nc, skip = FALSE), std)
    }
  }, error = function(e) {
    cat("  Model error:", conditionMessage(e), "\n")
    ins_results[[ins]] <<- list(n_tot = n_tot, n_ca = n_ca, n_nc = n_nc, skip = TRUE)
  })
  cat("\n")
}

# ==============================================================================
# 5-4  Stratified analysis — by race/ethnicity (3 groups)
# ==============================================================================
cat("--- Step 3: Stratified analysis by race/ethnicity ---\n\n")

race_levels  <- c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic")
race_results <- list()

for (rg in race_levels) {
  cat("Race/ethnicity:", rg, "\n")
  d_race <- subset(nhis_design, race_ethnicity == rg)
  n_tot  <- nrow(d_race$variables)
  n_ca   <- sum(d_race$variables$cancer_history == "Cancer history",    na.rm = TRUE)
  n_nc   <- sum(d_race$variables$cancer_history == "No cancer history", na.rm = TRUE)
  cat(sprintf("  n total = %d | Cancer = %d | No cancer = %d\n", n_tot, n_ca, n_nc))
  
  if (n_ca < NCHS_THRESHOLD || n_nc < NCHS_THRESHOLD) {
    cat(sprintf("  Suppressed (n < %d per group per NCHS guidelines)\n\n", NCHS_THRESHOLD))
    race_results[[rg]] <- list(n_tot = n_tot, n_ca = n_ca, n_nc = n_nc, skip = TRUE)
    next
  }
  
  tryCatch({
    m_race <- svyglm(
      glp1_use ~ cancer_history +
        age_years + sex +
        education + income_fpl_2cat + insurance_type + region +
        bmi_category + ascvd + ckd,
      design = d_race, family = quasibinomial()
    )
    std <- marg_std_subset(m_race, d_race, rg)
    if (is.null(std)) {
      race_results[[rg]] <- list(n_tot = n_tot, n_ca = n_ca, n_nc = n_nc, skip = TRUE)
    } else {
      cat(sprintf("  Cancer: %.1f%% (SE %.1f) | No cancer: %.1f%% (SE %.1f) | Gap: %.1f pp\n",
                  std$prev_ca, std$se_ca, std$prev_nc, std$se_nc, std$diff))
      race_results[[rg]] <- c(list(n_tot = n_tot, n_ca = n_ca, n_nc = n_nc, skip = FALSE), std)
    }
  }, error = function(e) {
    cat("  Model error:", conditionMessage(e), "\n")
    race_results[[rg]] <<- list(n_tot = n_tot, n_ca = n_ca, n_nc = n_nc, skip = TRUE)
  })
  cat("\n")
}

# ==============================================================================
# 5-5  Build Table 3
# ==============================================================================
cat("--- Step 4: Building Table 3 ---\n\n")

make_row <- function(label, res) {
  if (isTRUE(res$skip)) {
    data.frame(Stratum = label, n = sprintf("%d / %d", res$n_ca, res$n_nc),
               prev_ca = "\u2014", prev_nc = "\u2014", diff = "\u2014",
               stringsAsFactors = FALSE)
  } else {
    data.frame(Stratum = label, n = sprintf("%d / %d", res$n_ca, res$n_nc),
               prev_ca = fmt_prev(res$prev_ca, res$se_ca),
               prev_nc = fmt_prev(res$prev_nc, res$se_nc),
               diff    = fmt_diff(res$diff, res$ci_lo, res$ci_hi),
               stringsAsFactors = FALSE)
  }
}

col_names <- c("Stratum", "n (Cancer / No cancer)",
               "Cancer history,\n% (SE)", "No cancer history,\n% (SE)",
               "Absolute difference\n(Cancer \u2212 No cancer), pp (95% CI)")

t3a <- bind_rows(lapply(insurance_levels, function(x) make_row(x, ins_results[[x]])))
t3b <- bind_rows(lapply(race_levels,      function(x) make_row(x, race_results[[x]])))
colnames(t3a) <- colnames(t3b) <- col_names

# ==============================================================================
# 5-6  Flextables and export
# ==============================================================================
make_ft3 <- function(df, panel_label, p_int) {
  flextable(df) %>%
    theme_booktabs() %>% autofit() %>%
    align(align = "left",   j = 1,   part = "all") %>%
    align(align = "center", j = 2:5, part = "body") %>%
    align(align = "center", j = 2:5, part = "header") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    fontsize(size = 10, part = "body") %>%
    fontsize(size = 10, part = "header") %>%
    bold(part = "header") %>%
    add_header_lines(
      values = sprintf("%s (Model B standardized estimates; P for interaction = %s)",
                       panel_label, fmt_p_val(p_int))
    )
}

ft3a <- make_ft3(t3a, "A. Stratified by Insurance Type",  p_int_ins)
ft3b <- make_ft3(t3b, "B. Stratified by Race/Ethnicity", p_int_race)

footnote_t3 <- paste(
  "Survey-weighted estimates from Model B (adjusted for age, sex, education,",
  "income [<200% vs \u2265200% FPL], region, BMI, ASCVD, and CKD) using marginal",
  "standardization within each stratum.",
  "Insurance models additionally adjust for race/ethnicity;",
  "race/ethnicity models additionally adjust for insurance type.",
  "Uninsured respondents (n = 182) excluded from all analyses.",
  "'Other' race/ethnicity category excluded from Panel B due to heterogeneity",
  "and small cell sizes in the cancer survivor stratum.",
  sprintf("Strata with fewer than %d unweighted observations in either cancer history",
          NCHS_THRESHOLD),
  "group are suppressed (\u2014), consistent with NCHS analytic guidelines.",
  "Absolute difference 95% CI computed via delta method.",
  sprintf("P for interaction (cancer \u00d7 insurance type): %s.", fmt_p_val(p_int_ins)),
  sprintf("P for interaction (cancer \u00d7 race/ethnicity): %s.", fmt_p_val(p_int_race)),
  "SE = standard error; pp = percentage points; CI = confidence interval."
)

fname <- ifelse(include_main,
                "table3_effect_modification_MAIN.docx",
                "table3_effect_modification_SUPPLEMENT.docx")

doc3 <- read_docx() %>%
  body_add_par(
    "Table 3. Effect Modification of the Cancer\u2013GLP-1 RA Association by Insurance Type and Race/Ethnicity (NHIS 2024)",
    style = "heading 1"
  ) %>%
  body_add_par("") %>%
  body_add_flextable(ft3a) %>%
  body_add_par("") %>%
  body_add_flextable(ft3b) %>%
  body_add_par("") %>%
  body_add_par(footnote_t3, style = "Normal")

print(doc3, target = fname)

cat("\n=============================================================================\n")
cat("Block 05 complete:", fname, "\n")
cat("  Cancer x Insurance P =",      fmt_p_val(p_int_ins), "\n")
cat("  Cancer x Race/Ethnicity P =", fmt_p_val(p_int_race), "\n")
cat("=============================================================================\n")


# ==============================================================================
# NHIS 2024 Analysis Pipeline — Block 06: eTable 1 — Regression Odds Ratios
# ==============================================================================
# Input:  model_a, model_b (Block 04)
# Output: eTable1_regression_coefficients.docx
#
# Extracts and formats survey-weighted logistic regression odds ratios (ORs)
# and 95% CIs for all predictors from Models A and B.
#
# Income: 2-category primary specification (<200% vs >=200% FPL).
# BMI, ASCVD, CKD appear in Model B only.
# ==============================================================================

library(dplyr)
library(survey)
library(flextable)
library(officer)

if (!exists("model_a")) stop("model_a not found. Run Block 04 first.")
if (!exists("model_b")) stop("model_b not found. Run Block 04 first.")

cat("=============================================================================\n")
cat("Block 06: eTable 1 \u2014 Regression Odds Ratios\n")
cat("=============================================================================\n\n")

# ==============================================================================
# 6-0  Helpers
# ==============================================================================
extract_coef <- function(model, coef_name, model_label = "") {
  coefs <- coef(model); vcovs <- vcov(model)
  summ  <- summary(model)$coefficients
  if (!coef_name %in% names(coefs)) {
    cat(sprintf("  WARNING: '%s' not found in %s model\n", coef_name, model_label))
    return(list(or = NA, lo = NA, hi = NA, p = NA))
  }
  est <- coefs[coef_name]; se <- sqrt(vcovs[coef_name, coef_name])
  list(or = exp(est), lo = exp(est - 1.96*se), hi = exp(est + 1.96*se),
       p = summ[coef_name, "Pr(>|t|)"])
}

fmt_or <- function(res) {
  if (any(is.na(c(res$or, res$lo, res$hi)))) return("\u2014")
  sprintf("%.2f (%.2f, %.2f)", res$or, res$lo, res$hi)
}
fmt_p <- function(p) {
  if (is.null(p) || is.na(p)) return("\u2014")
  if (p < 0.001) return("<.001")
  sprintf("%.3f", p)
}

row_ref       <- function(var, ref_label) c(var, paste0(ref_label, " (Ref)"), "\u2014","\u2014","\u2014","\u2014")
row_coef      <- function(var_label, cat_label, ca, cb) c(var_label, cat_label, fmt_or(ca), fmt_p(ca$p), fmt_or(cb), fmt_p(cb$p))
row_coef_b    <- function(var_label, cat_label, cb)   c(var_label, cat_label, "\u2014","\u2014", fmt_or(cb), fmt_p(cb$p))

# ==============================================================================
# 6-1  Extract all coefficients
# ==============================================================================
cat("Extracting coefficients...\n\n")

results <- list(); i <- 0L
add <- function(row) { i <<- i + 1L; results[[i]] <<- row }

add(row_ref("Cancer history", "No cancer history"))
add(row_coef("", "Cancer history",
             extract_coef(model_a, "cancer_historyCancer history", "A"),
             extract_coef(model_b, "cancer_historyCancer history", "B")))

add(row_coef("Age, years", "Per 1-year increase",
             extract_coef(model_a, "age_years", "A"),
             extract_coef(model_b, "age_years", "B")))

add(row_ref("Sex", "Male"))
add(row_coef("", "Female",
             extract_coef(model_a, "sexFemale", "A"),
             extract_coef(model_b, "sexFemale", "B")))

add(row_ref("Race/ethnicity", "Non-Hispanic White"))
for (lv in c("Non-Hispanic Black", "Hispanic", "Other")) {
  add(row_coef("", lv,
               extract_coef(model_a, paste0("race_ethnicity", lv), "A"),
               extract_coef(model_b, paste0("race_ethnicity", lv), "B")))
}

add(row_ref("Education", "High school or less"))
add(row_coef("", "Some college",
             extract_coef(model_a, "educationSome college", "A"),
             extract_coef(model_b, "educationSome college", "B")))
add(row_coef("", "Bachelor's+",
             extract_coef(model_a, "educationBachelor's+", "A"),
             extract_coef(model_b, "educationBachelor's+", "B")))

add(row_ref("Income (% FPL)", "<200% FPL"))
add(row_coef("", "\u2265200% FPL",
             extract_coef(model_a, "income_fpl_2cat\u2265200% FPL", "A"),
             extract_coef(model_b, "income_fpl_2cat\u2265200% FPL", "B")))

add(row_ref("Insurance type", "Private"))
add(row_coef("", "Medicare",
             extract_coef(model_a, "insurance_typeMedicare", "A"),
             extract_coef(model_b, "insurance_typeMedicare", "B")))
add(row_coef("", "Medicaid",
             extract_coef(model_a, "insurance_typeMedicaid", "A"),
             extract_coef(model_b, "insurance_typeMedicaid", "B")))

add(row_ref("Region", "Northeast"))
for (lv in c("Midwest", "South", "West")) {
  add(row_coef("", lv,
               extract_coef(model_a, paste0("region", lv), "A"),
               extract_coef(model_b, paste0("region", lv), "B")))
}

add(row_ref("BMI (kg/m\u00b2)", "<30"))
add(row_coef_b("", "30\u201334.9",
               extract_coef(model_b, "bmi_category30\u201334.9", "B")))
add(row_coef_b("", "\u226535",
               extract_coef(model_b, "bmi_category\u226535", "B")))

add(row_coef_b("ASCVD", "Yes vs No",
               extract_coef(model_b, "ascvdYes", "B")))
add(row_coef_b("Chronic kidney disease", "Yes vs No",
               extract_coef(model_b, "ckdYes", "B")))

cat("Coefficients extracted:", length(results), "rows\n\n")

# ==============================================================================
# 6-2  Build data frame and flextable
# ==============================================================================
et1_df <- as.data.frame(do.call(rbind, results), stringsAsFactors = FALSE)
colnames(et1_df) <- c("Predictor","Category",
                      "Model A OR (95% CI)","Model A P value",
                      "Model B OR (95% CI)","Model B P value")

header_idx <- which(et1_df$Predictor != "" &
                      (grepl("\\(Ref\\)$", et1_df$Category) |
                         et1_df$Category == "Per 1-year increase"))

ft_et1 <- flextable(et1_df) %>%
  theme_booktabs() %>% autofit() %>%
  align(align = "left",   j = 1:2, part = "all") %>%
  align(align = "center", j = 3:6, part = "body") %>%
  align(align = "center", j = 3:6, part = "header") %>%
  font(fontname = "Times New Roman", part = "all") %>%
  fontsize(size = 10, part = "body") %>%
  fontsize(size = 10, part = "header") %>%
  bold(part = "header") %>%
  bold(i = header_idx, j = 1) %>%
  merge_v(j = 1) %>%
  valign(j = 1, valign = "top")

# ==============================================================================
# 6-3  Export to Word
# ==============================================================================
n_uninsured <- sum(is.na(analytic_df$insurance_type))

footnote_et1 <- paste(
  "Odds ratios (95% confidence intervals) from survey-weighted logistic regression",
  "accounting for the NHIS 2024 complex sampling design (PSUs, strata, annual weights).",
  "Model A: adjusted for cancer history, age (continuous), sex, race/ethnicity (4 categories),",
  "education (3 categories), income (<200% vs \u2265200% FPL), insurance type",
  "(Private/Medicare/Medicaid), and region.",
  paste0("Uninsured respondents (n = ", n_uninsured, ") excluded from all analyses."),
  "Model B: Model A plus BMI category (3 groups: <30, 30\u201334.9, \u226535 kg/m\u00b2),",
  "atherosclerotic cardiovascular disease (CHD/MI or stroke), and chronic kidney disease.",
  "Income shown as 2-category (primary specification); see eTable 2 (S1) for",
  "3-category income sensitivity analysis.",
  "'Other' race/ethnicity includes non-Hispanic Asian, non-Hispanic Pacific Islander,",
  "American Indian or Alaska Native, other single race, and multiracial respondents;",
  "collapsed due to sparse cells in the cancer survivor stratum.",
  "BMI collapsed to 3 categories due to sparse cells.",
  "\u2014 = not applicable (variable not included in this model).",
  "OR = odds ratio; CI = confidence interval; FPL = federal poverty level;",
  "ASCVD = atherosclerotic cardiovascular disease; CKD = chronic kidney disease."
)

doc_et1 <- read_docx() %>%
  body_add_par(
    "eTable 1. Survey-Weighted Logistic Regression Odds Ratios for Models A and B: GLP-1 RA Use Among U.S. Adults With Diabetes (NHIS 2024)",
    style = "heading 1"
  ) %>%
  body_add_par("") %>%
  body_add_flextable(ft_et1) %>%
  body_add_par("") %>%
  body_add_par(footnote_et1, style = "Normal")

print(doc_et1, target = "eTable1_regression_coefficients.docx")

cat("\n=============================================================================\n")
cat("Block 06 complete \u2014 eTable1_regression_coefficients.docx saved\n")
cat("=============================================================================\n")

# ==============================================================================
# NHIS 2024 Analysis Pipeline — Block 07: eTable 2 — Sensitivity Analyses
# ==============================================================================
# Input:  analytic_df (Block 01), nhis_design (Block 02), model_b (Block 04)
# Output: eTable2_sensitivity_analyses.docx
#
# Sensitivity analyses:
#   Primary — Model B reference
#   S1      — Income 3-category (<200%, 200-399%, >=400% FPL)
#   S2      — Model B + hypertension + hyperlipidemia
#   S3      — Restrict to >=1 outpatient visit in past 12 months
#   S4      — Insurance 2-category (Private vs Public)
#   S5      — Subgroup stability (age group, sex, BMI; Model A specification)
#   S6      — Restrict to adults with type 2 diabetes only
#
# Cell suppression threshold: n >= 30 per cancer history group (NCHS guidelines)
# ==============================================================================

library(dplyr)
library(survey)
library(flextable)
library(officer)

if (!exists("nhis_design")) stop("nhis_design not found. Run Block 02 first.")
if (!exists("model_b"))     stop("model_b not found. Run Block 04 first.")
if (!exists("analytic_df")) stop("analytic_df not found. Run Block 01 first.")

NCHS_THRESHOLD <- 30

cat("=============================================================================\n")
cat("Block 07: eTable 2 -- Sensitivity Analyses (S1-S6)\n")
cat(sprintf("NCHS cell suppression threshold: n >= %d\n", NCHS_THRESHOLD))
cat("=============================================================================\n\n")

# ==============================================================================
# 7-0  Formatting helpers
# ==============================================================================
fmt_gap    <- function(g, lo, hi) { if (any(is.na(c(g,lo,hi)))) return("\u2014"); sprintf("%.1f (%.1f, %.1f)", g, lo, hi) }
fmt_or     <- function(or, lo, hi){ if (any(is.na(c(or,lo,hi)))) return("\u2014"); sprintf("%.2f (%.2f, %.2f)", or, lo, hi) }
format_n   <- function(n) { if (is.na(n)) return("NA"); formatC(as.integer(n), format="d", big.mark=",") }

# ==============================================================================
# 7-1  Sensitivity extraction helper
# ==============================================================================
extract_sens <- function(model, design, label) {
  cat("  Extracting:", label, "\n")
  tryCatch({
    mv    <- intersect(all.vars(formula(model)), names(design$variables))
    design$variables$cc_sens <- complete.cases(design$variables[, mv])
    design_cc <- subset(design, cc_sens)
    df_cc     <- design_cc$variables
    n_ca      <- sum(df_cc$cancer_history == "Cancer history",    na.rm = TRUE)
    n_nc      <- sum(df_cc$cancer_history == "No cancer history", na.rm = TRUE)
    n_tot     <- nrow(df_cc)
    
    design_cc$variables$pred_ca <- as.numeric(predict(model,
                                                      newdata = df_cc %>% mutate(cancer_history = factor("Cancer history",    levels = levels(analytic_df$cancer_history))),
                                                      type = "response"))
    design_cc$variables$pred_nc <- as.numeric(predict(model,
                                                      newdata = df_cc %>% mutate(cancer_history = factor("No cancer history", levels = levels(analytic_df$cancer_history))),
                                                      type = "response"))
    
    obj_ca <- svymean(~pred_ca, design_cc, na.rm = TRUE)
    obj_nc <- svymean(~pred_nc, design_cc, na.rm = TRUE)
    p_ca   <- 100 * as.numeric(coef(obj_ca)); se_ca <- 100 * as.numeric(SE(obj_ca))
    p_nc   <- 100 * as.numeric(coef(obj_nc)); se_nc <- 100 * as.numeric(SE(obj_nc))
    gap    <- p_ca - p_nc; se_gap <- sqrt(se_ca^2 + se_nc^2)
    
    coef_ca <- coef(model)["cancer_historyCancer history"]
    se_or   <- sqrt(vcov(model)["cancer_historyCancer history","cancer_historyCancer history"])
    
    list(gap=gap, gap_lo=gap-1.96*se_gap, gap_hi=gap+1.96*se_gap,
         or=exp(coef_ca), or_lo=exp(coef_ca-1.96*se_or), or_hi=exp(coef_ca+1.96*se_or),
         n_tot=n_tot, n_ca=n_ca, n_nc=n_nc, ok=TRUE)
  }, error = function(e) {
    cat("  ERROR:", conditionMessage(e), "\n")
    list(gap=NA, gap_lo=NA, gap_hi=NA, or=NA, or_lo=NA, or_hi=NA,
         n_tot=NA, n_ca=NA, n_nc=NA, ok=FALSE)
  })
}

# ==============================================================================
# 7-2  Primary (Model B)
# ==============================================================================
cat("--- Primary: Model B ---\n")
prim <- extract_sens(model_b, nhis_design, "Primary Model B")
cat(sprintf("  Gap = %.1f pp (%.1f, %.1f) | aOR = %.2f (%.2f, %.2f) | n = %s\n\n",
            prim$gap, prim$gap_lo, prim$gap_hi, prim$or, prim$or_lo, prim$or_hi, format_n(prim$n_tot)))

# ==============================================================================
# 7-3  S1: Income 3-category
# ==============================================================================
cat("--- S1: Income 3-category ---\n")
nhis_design$variables$income_fpl_3cat <- factor(
  case_when(
    nhis_design$variables$income_fpl_4cat %in% c("<100% FPL","100\u2013199% FPL") ~ "<200% FPL",
    nhis_design$variables$income_fpl_4cat == "200\u2013399% FPL" ~ "200\u2013399% FPL",
    nhis_design$variables$income_fpl_4cat == "\u2265400% FPL"    ~ "\u2265400% FPL",
    TRUE ~ NA_character_
  ), levels = c("<200% FPL","200\u2013399% FPL","\u2265400% FPL")
)
model_s1 <- svyglm(
  glp1_use ~ cancer_history + age_years + sex + race_ethnicity +
    education + income_fpl_3cat + insurance_type + region + bmi_category + ascvd + ckd,
  design = nhis_design, family = quasibinomial()
)
s1 <- extract_sens(model_s1, nhis_design, "S1")
cat(sprintf("  Gap = %.1f pp (%.1f, %.1f) | aOR = %.2f (%.2f, %.2f)\n\n",
            s1$gap, s1$gap_lo, s1$gap_hi, s1$or, s1$or_lo, s1$or_hi))

# ==============================================================================
# 7-4  S2: Model B + hypertension + hyperlipidemia
# ==============================================================================
cat("--- S2: Model B + hypertension + hyperlipidemia ---\n")
model_s2 <- svyglm(
  glp1_use ~ cancer_history + age_years + sex + race_ethnicity +
    education + income_fpl_2cat + insurance_type + region +
    bmi_category + ascvd + ckd + hypertension + hyperlipidemia,
  design = nhis_design, family = quasibinomial()
)
s2 <- extract_sens(model_s2, nhis_design, "S2")
cat(sprintf("  Gap = %.1f pp (%.1f, %.1f) | aOR = %.2f (%.2f, %.2f)\n\n",
            s2$gap, s2$gap_lo, s2$gap_hi, s2$or, s2$or_lo, s2$or_hi))

# ==============================================================================
# 7-5  S3: Restrict to >=1 outpatient visit
# ==============================================================================
cat("--- S3: Restrict to outpatient visit past 12 months ---\n")
design_s3 <- subset(nhis_design, outpatient_12m == "Yes")
cat(sprintf("  n before: %s | after: %s (%.1f%%)\n",
            format_n(nrow(nhis_design$variables)), format_n(nrow(design_s3$variables)),
            100 * nrow(design_s3$variables) / nrow(nhis_design$variables)))
model_s3 <- svyglm(
  glp1_use ~ cancer_history + age_years + sex + race_ethnicity +
    education + income_fpl_2cat + insurance_type + region + bmi_category + ascvd + ckd,
  design = design_s3, family = quasibinomial()
)
s3 <- extract_sens(model_s3, design_s3, "S3")
cat(sprintf("  Gap = %.1f pp (%.1f, %.1f) | aOR = %.2f (%.2f, %.2f)\n\n",
            s3$gap, s3$gap_lo, s3$gap_hi, s3$or, s3$or_lo, s3$or_hi))

# ==============================================================================
# 7-6  S4: Insurance 2-category (Private vs Public)
# ==============================================================================
cat("--- S4: Insurance 2-category ---\n")
nhis_design$variables$insurance_2cat <- factor(
  case_when(
    nhis_design$variables$insurance_type == "Private"                    ~ "Private",
    nhis_design$variables$insurance_type %in% c("Medicare","Medicaid")  ~ "Public",
    TRUE ~ NA_character_
  ), levels = c("Private","Public")
)
model_s4 <- svyglm(
  glp1_use ~ cancer_history + age_years + sex + race_ethnicity +
    education + income_fpl_2cat + insurance_2cat + region + bmi_category + ascvd + ckd,
  design = nhis_design, family = quasibinomial()
)
s4 <- extract_sens(model_s4, nhis_design, "S4")
cat(sprintf("  Gap = %.1f pp (%.1f, %.1f) | aOR = %.2f (%.2f, %.2f)\n\n",
            s4$gap, s4$gap_lo, s4$gap_hi, s4$or, s4$or_lo, s4$or_hi))

# ==============================================================================
# 7-7  S5: Subgroup stability (Model A specification; NCHS threshold n >= 30)
# ==============================================================================
cat("--- S5: Subgroup stability ---\n\n")

subgroups <- list(
  list(var="age_group",    levels=c("18\u201344","45\u201364","65+")),
  list(var="sex",          levels=c("Male","Female")),
  list(var="bmi_category", levels=c("<30","30\u201334.9","\u226535"))
)

s5_rows <- list()

for (sg in subgroups) {
  cat("  Subgroup:", sg$var, "\n")
  for (lv in sg$levels) {
    d_sub <- tryCatch(subset(nhis_design, get(sg$var) == lv), error = function(e) NULL)
    if (is.null(d_sub)) {
      s5_rows[[length(s5_rows)+1]] <- list(sg$var, lv, "\u2014", "n=NA", "Subset error")
      next
    }
    n_tot <- nrow(d_sub$variables)
    n_ca  <- sum(d_sub$variables$cancer_history == "Cancer history",    na.rm=TRUE)
    n_nc  <- sum(d_sub$variables$cancer_history == "No cancer history", na.rm=TRUE)
    cat(sprintf("    %s: n = %s (Cancer = %d, No cancer = %d) ", lv, format_n(n_tot), n_ca, n_nc))
    
    if (n_ca < NCHS_THRESHOLD || n_nc < NCHS_THRESHOLD) {
      cat(sprintf("-- suppressed (n < %d per NCHS guidelines)\n", NCHS_THRESHOLD))
      s5_rows[[length(s5_rows)+1]] <- list(sg$var, lv, "\u2014",
                                           sprintf("n = %s", format_n(n_tot)), "Insufficient n")
      next
    }
    tryCatch({
      m_sub <- svyglm(
        glp1_use ~ cancer_history + age_years + sex + race_ethnicity +
          education + income_fpl_2cat + insurance_type + region,
        design = d_sub, family = quasibinomial()
      )
      coef_ca <- coef(m_sub)["cancer_historyCancer history"]
      se_ca   <- sqrt(vcov(m_sub)["cancer_historyCancer history","cancer_historyCancer history"])
      or <- exp(coef_ca); lo <- exp(coef_ca-1.96*se_ca); hi <- exp(coef_ca+1.96*se_ca)
      cat(sprintf("-- aOR = %.2f (%.2f, %.2f)\n", or, lo, hi))
      s5_rows[[length(s5_rows)+1]] <- list(sg$var, lv,
                                           sprintf("%.2f (%.2f, %.2f)", or, lo, hi),
                                           sprintf("n = %s", format_n(n_tot)), "")
    }, error = function(e) {
      cat("-- model error:", conditionMessage(e), "\n")
      s5_rows[[length(s5_rows)+1]] <- list(sg$var, lv, "\u2014",
                                           sprintf("n = %s", format_n(n_tot)), "Model error")
    })
  }
  cat("\n")
}

# ==============================================================================
# 7-8  S6: Restrict to type 2 diabetes only
# ==============================================================================
cat("--- S6: Restrict to adults with type 2 diabetes only ---\n")
nhis_design$variables$type2_flag <- analytic_df$DIBTYPE_A == 2
design_s6 <- subset(nhis_design, type2_flag)

pct_ca_t1 <- 100 * sum(analytic_df$DIBTYPE_A==1 & analytic_df$cancer_history=="Cancer history", na.rm=TRUE) /
  sum(analytic_df$cancer_history=="Cancer history", na.rm=TRUE)
pct_nc_t1 <- 100 * sum(analytic_df$DIBTYPE_A==1 & analytic_df$cancer_history=="No cancer history", na.rm=TRUE) /
  sum(analytic_df$cancer_history=="No cancer history", na.rm=TRUE)
cat(sprintf("  Type 1 proportion: Cancer = %.1f%% | No cancer = %.1f%%\n", pct_ca_t1, pct_nc_t1))
cat(sprintf("  Type 2 subset n: %s\n", format_n(nrow(design_s6$variables))))

model_s6 <- svyglm(
  glp1_use ~ cancer_history + age_years + sex + race_ethnicity +
    education + income_fpl_2cat + insurance_type + region + bmi_category + ascvd + ckd,
  design = design_s6, family = quasibinomial()
)
s6 <- extract_sens(model_s6, design_s6, "S6: Type 2 diabetes only")
if (!is.na(s6$or)) {
  cat(sprintf("  Gap = %.1f pp (%.1f, %.1f) | aOR = %.2f (%.2f, %.2f) | n = %s\n",
              s6$gap, s6$gap_lo, s6$gap_hi, s6$or, s6$or_lo, s6$or_hi, format_n(s6$n_tot)))
  cat(sprintf("  %s\n\n", ifelse(s6$or_lo <= 1.0 & s6$or_hi >= 1.0,
                                 "NOT statistically significant (95% CI includes 1.0)", "statistically significant")))
}

# ==============================================================================
# 7-9  Build eTable 2
# ==============================================================================
cat("--- Building eTable 2 ---\n\n")

s6_sig  <- ifelse(!is.na(s6$or_lo) && s6$or_lo <= 1.0,
                  "directionally consistent; not statistically significant (95% CI includes 1.0)",
                  "directionally consistent; statistically significant")

et2_df <- data.frame(
  Analysis = c("Primary","S1","S2","S3","S4","S5","S6"),
  Modification = c(
    "Model B (reference)",
    "Income 3-category (<200%, 200\u2013399%, \u2265400% FPL)",
    "Model B + hypertension + hyperlipidemia",
    "Restrict to \u22651 outpatient visit in past 12 months",
    "Insurance 2-category (Private vs Public [Medicare+Medicaid])",
    "Subgroup stability (see panel below)",
    "Restrict to adults with type 2 diabetes only"
  ),
  Gap = c(fmt_gap(prim$gap,prim$gap_lo,prim$gap_hi),
          fmt_gap(s1$gap,s1$gap_lo,s1$gap_hi),
          fmt_gap(s2$gap,s2$gap_lo,s2$gap_hi),
          fmt_gap(s3$gap,s3$gap_lo,s3$gap_hi),
          fmt_gap(s4$gap,s4$gap_lo,s4$gap_hi),
          "\u2014",
          fmt_gap(s6$gap,s6$gap_lo,s6$gap_hi)),
  OR  = c(fmt_or(prim$or,prim$or_lo,prim$or_hi),
          fmt_or(s1$or,s1$or_lo,s1$or_hi),
          fmt_or(s2$or,s2$or_lo,s2$or_hi),
          fmt_or(s3$or,s3$or_lo,s3$or_hi),
          fmt_or(s4$or,s4$or_lo,s4$or_hi),
          "\u2014",
          fmt_or(s6$or,s6$or_lo,s6$or_hi)),
  Notes = c(
    sprintf("n = %s (reference)", format_n(prim$n_tot)),
    sprintf("n = %s", format_n(s1$n_tot)),
    sprintf("n = %s", format_n(s2$n_tot)),
    sprintf("n = %s (%.1f%% of primary)", format_n(s3$n_tot), 100*s3$n_tot/prim$n_tot),
    sprintf("n = %s", format_n(s4$n_tot)),
    "See subgroup panel",
    sprintf("n = %s; %s", format_n(s6$n_tot), s6_sig)
  ),
  stringsAsFactors = FALSE
)
colnames(et2_df) <- c("Analysis","Modification","Standardized gap, pp\n(95% CI)","aOR\n(95% CI)","Notes")

ft_et2 <- flextable(et2_df) %>%
  theme_booktabs() %>% autofit() %>%
  align(align="left",   j=1:2, part="all") %>%
  align(align="center", j=3:5, part="body") %>%
  align(align="center", j=3:5, part="header") %>%
  font(fontname="Times New Roman", part="all") %>%
  fontsize(size=10, part="body") %>%
  fontsize(size=10, part="header") %>%
  bold(part="header") %>% bold(i=1, j=1) %>%
  bg(i=7, bg="#F0F4FB")

# S5 subgroup panel
if (length(s5_rows) == 0) {
  s5_df <- data.frame(Subgroup="No subgroup results", Level="\u2014",
                      `aOR (95% CI)`="\u2014", n="\u2014", Note="See console output",
                      stringsAsFactors=FALSE, check.names=FALSE)
} else {
  s5_df <- as.data.frame(do.call(rbind, lapply(s5_rows, function(x) as.character(unlist(x)))),
                         stringsAsFactors=FALSE)
  colnames(s5_df) <- c("Subgroup","Level","aOR (95% CI)","n","Note")
}

ft_s5 <- flextable(s5_df) %>%
  theme_booktabs() %>% autofit() %>%
  align(align="left",   j=1:2, part="all") %>%
  align(align="center", j=3:5, part="body") %>%
  align(align="center", j=3:5, part="header") %>%
  font(fontname="Times New Roman", part="all") %>%
  fontsize(size=10, part="body") %>%
  fontsize(size=10, part="header") %>%
  bold(part="header") %>%
  merge_v(j=1) %>% valign(j=1, valign="top") %>%
  add_header_lines("S5. Subgroup Stability: Cancer\u2013GLP-1 RA aOR Within Subgroups")

# Footnote
footnote_et2 <- paste0(
  "All sensitivity analyses use survey-weighted logistic regression with marginal standardization, ",
  "accounting for the NHIS 2024 complex sampling design (PSUs, strata, annual weights). ",
  "Primary (Model B): adjusted for cancer history, age (continuous), sex, race/ethnicity ",
  "(4 categories), education (3 categories), income (<200% vs \u2265200% FPL), insurance type ",
  "(Private/Medicare/Medicaid), region, BMI (3 categories: <30, 30\u201334.9, \u226535\u00a0kg/m\u00b2), ASCVD, and CKD. ",
  "Uninsured respondents (n\u00a0=\u00a0182) excluded from all analyses. ",
  "S1: income recategorized to 3 groups (<200%, 200\u2013399%, \u2265400% FPL). ",
  "S2: Model B + hypertension + hyperlipidemia. ",
  "S3: restricted to respondents with \u22651 outpatient visit in the past 12 months. ",
  "S4: insurance collapsed to 2 categories (Private vs Public = Medicare + Medicaid). ",
  "S5: subgroup-specific aORs from Model A specification (without BMI, ASCVD, and CKD) due to ",
  "smaller stratum sample sizes; strata with fewer than 30 unweighted observations per cancer ",
  "history group are suppressed, consistent with NCHS analytic guidelines for NHIS public-use data. ",
  "S6: restricted to adults with type 2 diabetes only, to assess robustness to diabetes type ",
  "classification. The primary analysis follows NCHS convention by including all adults with ",
  "diagnosed diabetes; type 1 diabetes accounted for 8.9% of cancer survivors and 8.8% of those ",
  "without cancer history, confirming comparable distributions across groups. ",
  "Standardized gap = percentage-point difference (Cancer \u2212 No cancer); 95% CI by delta method. ",
  "aOR = adjusted odds ratio; pp = percentage points; FPL = federal poverty level; ",
  "ASCVD = atherosclerotic cardiovascular disease; CKD = chronic kidney disease."
)

doc_et2 <- read_docx() %>%
  body_add_par("eTable 2. Sensitivity Analyses for the Association Between Cancer History and GLP-1 RA Use (NHIS 2024)",
               style="heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft_et2) %>%
  body_add_par("") %>%
  body_add_flextable(ft_s5) %>%
  body_add_par("") %>%
  body_add_par(footnote_et2, style="Normal")

print(doc_et2, target="eTable2_sensitivity_analyses.docx")

cat("\n=============================================================================\n")
cat("Block 07 RESULTS SUMMARY\n")
cat("=============================================================================\n")
for (nm in c("Primary","S1","S2","S3","S4","S6")) {
  res <- list(Primary=prim, S1=s1, S2=s2, S3=s3, S4=s4, S6=s6)[[nm]]
  if (!is.na(res$or)) {
    sig <- ifelse(res$or_lo > 1.0 | res$or_hi < 1.0, "significant", "NOT significant (CI includes 1.0)")
    cat(sprintf("  %-8s aOR = %.2f (%.2f, %.2f)  Gap = %.1f pp  n = %s  [%s]\n",
                nm, res$or, res$or_lo, res$or_hi, res$gap, format_n(res$n_tot), sig))
  }
}
cat("  S5       See subgroup panel\n")
cat("=============================================================================\n")
cat("Block 07 complete -- eTable2_sensitivity_analyses.docx saved\n")
cat("=============================================================================\n")


# ==============================================================================
# NHIS 2024 Analysis Pipeline — Block 08: eTable 3 — Exploratory Attenuation
# ==============================================================================
# Input:  analytic_df (Block 01), nhis_design (Block 02), model_b (Block 04)
# Output: eTable3_model_c_attenuation.docx
#         model_c  (used in Block 10)
#
# Fits Model C = Model B + cost_barrier + spd_k6
# Assesses whether these pathway variables attenuate the cancer-GLP-1 RA
# association. Attenuation = % change in cancer history log-OR from B to C.
#
# cost_barrier: any delayed or unmet care due to cost in past 12 months
# spd_k6:       Kessler-6 serious psychological distress indicator
#
# Model C is exploratory attenuation only — NOT causal mediation.
# Cross-sectional design precludes establishing temporality.
# ==============================================================================

library(dplyr)
library(survey)
library(flextable)
library(officer)

if (!exists("nhis_design")) stop("nhis_design not found. Run Block 02 first.")
if (!exists("model_b"))     stop("model_b not found. Run Block 04 first.")

cat("=============================================================================\n")
cat("Block 08: eTable 3 \u2014 Exploratory Attenuation Model (Model C)\n")
cat("=============================================================================\n\n")

# Confirm required variables exist in analytic_df
for (v in c("cost_barrier","spd_k6")) {
  if (!v %in% names(analytic_df))
    stop(sprintf("'%s' not found in analytic_df. Ensure Block 01 has been run.", v))
}

cat("cost_barrier distribution:\n"); print(table(analytic_df$cost_barrier, useNA="ifany"))
cat("spd_k6 distribution:\n");       print(table(analytic_df$spd_k6,       useNA="ifany"))
cat("\n")

# Inject variables into design (created after nhis_design was built)
nhis_design$variables$cost_barrier <- analytic_df$cost_barrier
nhis_design$variables$spd_k6       <- analytic_df$spd_k6

# Fit Model C
cat("Fitting Model C...\n")
model_c <- svyglm(
  glp1_use ~ cancer_history + age_years + sex + race_ethnicity +
    education + income_fpl_2cat + insurance_type + region +
    bmi_category + ascvd + ckd + cost_barrier + spd_k6,
  design = nhis_design, family = quasibinomial()
)
cat("Model C fitted\n\n")

# Extract ORs
extract_coef <- function(model, coef_name) {
  coefs <- coef(model); vcovs <- vcov(model); summ <- summary(model)$coefficients
  if (!coef_name %in% names(coefs)) {
    cat(sprintf("  WARNING: '%s' not found\n", coef_name))
    return(list(or=NA, lo=NA, hi=NA, p=NA))
  }
  est <- coefs[coef_name]; se <- sqrt(vcovs[coef_name, coef_name])
  list(or=exp(est), lo=exp(est-1.96*se), hi=exp(est+1.96*se), p=summ[coef_name,"Pr(>|t|)"])
}
fmt_or <- function(res) { if (any(is.na(c(res$or,res$lo,res$hi)))) return("\u2014"); sprintf("%.2f (%.2f, %.2f)", res$or, res$lo, res$hi) }
fmt_p  <- function(p)   { if (is.null(p)||is.na(p)) return("\u2014"); if (p<0.001) return("<.001"); sprintf("%.3f", p) }

ca_b    <- extract_coef(model_b, "cancer_historyCancer history")
ca_c    <- extract_coef(model_c, "cancer_historyCancer history")
cost_c  <- extract_coef(model_c, "cost_barrierYes")
spd_c   <- extract_coef(model_c, "spd_k6Yes")

pct_attenuation <- 100 * (log(ca_c$or) - log(ca_b$or)) / abs(log(ca_b$or))

cat("--- Attenuation summary ---\n")
cat(sprintf("  Model B cancer aOR: %.2f (%.2f, %.2f)  P = %s\n", ca_b$or, ca_b$lo, ca_b$hi, fmt_p(ca_b$p)))
cat(sprintf("  Model C cancer aOR: %.2f (%.2f, %.2f)  P = %s\n", ca_c$or, ca_c$lo, ca_c$hi, fmt_p(ca_c$p)))
cat(sprintf("  Attenuation:       %.1f%% (log-OR scale)\n", pct_attenuation))
cat(sprintf("  Cost barrier OR:   %.2f (%.2f, %.2f)  P = %s\n", cost_c$or, cost_c$lo, cost_c$hi, fmt_p(cost_c$p)))
cat(sprintf("  SPD (K6) OR:       %.2f (%.2f, %.2f)  P = %s\n\n", spd_c$or, spd_c$lo, spd_c$hi, fmt_p(spd_c$p)))

# Complete-case n
vars_b <- c("glp1_use","cancer_history","age_years","sex","race_ethnicity",
            "education","income_fpl_2cat","insurance_type","region","bmi_category","ascvd","ckd")
vars_c <- c(vars_b, "cost_barrier","spd_k6")
n_b <- sum(complete.cases(analytic_df[, vars_b]))
n_c <- sum(complete.cases(analytic_df[, vars_c]))
cat(sprintf("Complete cases: Model B n = %d | Model C n = %d\n\n", n_b, n_c))

# Build eTable 3
et3_df <- data.frame(
  Predictor = c("Cancer history","Cost barrier",
                "Serious psychological distress (K6)",
                "Other covariates","Sample size (complete case)"),
  Category  = c("Cancer vs No cancer history",
                "Any unmet or delayed care due to cost (Yes vs No)",
                "Kessler-6 score indicating serious psychological distress (Yes vs No)",
                "Age, sex, race/ethnicity, education, income, insurance, region, BMI, ASCVD, CKD",""),
  OR_B  = c(fmt_or(ca_b),"\u2014","\u2014","Included (see eTable 1)",sprintf("n = %s", formatC(n_b, format="d", big.mark=","))),
  p_B   = c(fmt_p(ca_b$p),"\u2014","\u2014","\u2014","\u2014"),
  OR_C  = c(fmt_or(ca_c),fmt_or(cost_c),fmt_or(spd_c),"Included (see eTable 1)",sprintf("n = %s", formatC(n_c, format="d", big.mark=","))),
  p_C   = c(fmt_p(ca_c$p),fmt_p(cost_c$p),fmt_p(spd_c$p),"\u2014","\u2014"),
  stringsAsFactors = FALSE
)
colnames(et3_df) <- c("Predictor","Category","Model B OR (95% CI)","Model B P value","Model C OR (95% CI)","Model C P value")

ft_et3 <- flextable(et3_df) %>%
  theme_booktabs() %>% autofit() %>%
  align(align="left",   j=1:2, part="all") %>%
  align(align="center", j=3:6, part="body") %>%
  align(align="center", j=3:6, part="header") %>%
  font(fontname="Times New Roman", part="all") %>%
  fontsize(size=10, part="body") %>%
  fontsize(size=10, part="header") %>%
  bold(part="header") %>%
  width(j=2, width=2.5)

attenuation_interp <- case_when(
  abs(pct_attenuation) < 5  ~ "minimal (\u22645%): cost barriers and psychological distress do not meaningfully explain the cancer\u2013GLP-1 RA association",
  abs(pct_attenuation) < 15 ~ "modest (5\u201315%): cost barriers and psychological distress partially explain the association",
  TRUE                      ~ "substantial (>15%): cost barriers and/or psychological distress are major contributors"
)

footnote_et3 <- paste(
  "Model C adds cost barrier and serious psychological distress to Model B as exploratory pathway variables.",
  "Model B: adjusted for age, sex, race/ethnicity, education, income (<200% vs \u2265200% FPL),",
  "insurance type (Private/Medicare/Medicaid), region, BMI, ASCVD, and CKD.",
  "Model C: Model B plus cost barrier and serious psychological distress.",
  "Cost barrier was defined as any delayed or unmet care due to cost in the past 12 months.",
  "Serious psychological distress was defined using the Kessler-6 scale.",
  sprintf("Attenuation of cancer history OR from Model B to Model C: %.1f%% (log-OR scale).", pct_attenuation),
  sprintf("Interpretation: attenuation is %s.", attenuation_interp),
  "Model C is presented as exploratory attenuation, not causal mediation.",
  "The cross-sectional design precludes establishing temporality or ruling out confounding;",
  "cost barrier and psychological distress may themselves be influenced by cancer history,",
  "and attenuation may therefore reflect partial overadjustment rather than mediation.",
  "OR = odds ratio; CI = confidence interval; SPD = serious psychological distress;",
  "ASCVD = atherosclerotic cardiovascular disease; CKD = chronic kidney disease."
)

doc_et3 <- read_docx() %>%
  body_add_par("eTable 3. Exploratory Attenuation Model (Model C): Cancer History, Cost Barrier, and Serious Psychological Distress as Predictors of GLP-1 RA Use (NHIS 2024)",
               style="heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft_et3) %>%
  body_add_par("") %>%
  body_add_par(footnote_et3, style="Normal")

print(doc_et3, target="eTable3_model_c_attenuation.docx")

cat("\n=============================================================================\n")
cat("Block 08 complete \u2014 eTable3_model_c_attenuation.docx saved\n")
cat(sprintf("  Attenuation: %.1f%% \u2014 %s\n", pct_attenuation,
            ifelse(abs(pct_attenuation)<5,"MINIMAL",ifelse(abs(pct_attenuation)<15,"MODEST","SUBSTANTIAL"))))
cat("=============================================================================\n")

# ==============================================================================
# NHIS 2024 Analysis Pipeline — Block 09: Figure 1
# ==============================================================================
# Input:  analytic_df (Block 01), nhis_design (Block 02), model_b (Block 04)
#         p_int_ins, p_int_race (Block 05)
# Output: Figure1_adjusted_probabilities.png / .pdf
#
# Panel A: Overall Model B standardized predicted probabilities.
# Panel B: By insurance type (shown if P for interaction < 0.10 from Block 05).
# Panel C: By race/ethnicity (shown if P for interaction < 0.10 from Block 05).
#
# All estimates use marginal standardization with survey-weighted svymean.
# Error bars represent 95% confidence intervals.
# ==============================================================================

library(dplyr)
library(survey)
library(ggplot2)

if (!exists("nhis_design")) stop("nhis_design not found. Run Block 02 first.")
if (!exists("model_b"))     stop("model_b not found. Run Block 04 first.")
if (!exists("p_int_ins"))   { cat("WARNING: p_int_ins not found; defaulting to 1.0\n");  p_int_ins  <- 1.0 }
if (!exists("p_int_race"))  { cat("WARNING: p_int_race not found; defaulting to 1.0\n"); p_int_race <- 1.0 }

THRESHOLD       <- 0.10
include_panel_b <- p_int_ins  < THRESHOLD
include_panel_c <- p_int_race < THRESHOLD

cat("=============================================================================\n")
cat("Block 09: Figure 1 \u2014 Adjusted Predicted Probabilities\n")
cat("=============================================================================\n\n")
cat(sprintf("  Panel B (Insurance):  P = %.4f  -> %s\n", p_int_ins,  ifelse(include_panel_b,"INCLUDE","OMIT")))
cat(sprintf("  Panel C (Race/Eth):   P = %.4f  -> %s\n\n", p_int_race, ifelse(include_panel_c,"INCLUDE","OMIT")))

# ==============================================================================
# 9-0  Marginal standardization helpers
# ==============================================================================
vars_b <- c("glp1_use","cancer_history","age_years","sex","race_ethnicity",
            "education","income_fpl_2cat","insurance_type","region","bmi_category","ascvd","ckd")

calc_overall <- function(model, design) {
  design$variables$cc_fig <- complete.cases(design$variables[, vars_b])
  d_cc <- subset(design, cc_fig); df_cc <- d_cc$variables
  
  d_cc$variables$pred_ca <- as.numeric(predict(model,
                                               newdata = df_cc %>% mutate(cancer_history = factor("Cancer history",    levels=levels(analytic_df$cancer_history))),
                                               type="response"))
  d_cc$variables$pred_nc <- as.numeric(predict(model,
                                               newdata = df_cc %>% mutate(cancer_history = factor("No cancer history", levels=levels(analytic_df$cancer_history))),
                                               type="response"))
  
  obj_ca <- svymean(~pred_ca, d_cc, na.rm=TRUE); obj_nc <- svymean(~pred_nc, d_cc, na.rm=TRUE)
  p_ca <- 100*as.numeric(coef(obj_ca)); se_ca <- 100*as.numeric(SE(obj_ca))
  p_nc <- 100*as.numeric(coef(obj_nc)); se_nc <- 100*as.numeric(SE(obj_nc))
  diff <- p_ca - p_nc; se_diff <- sqrt(se_ca^2 + se_nc^2)
  list(ca_est=p_ca, ca_lo=p_ca-1.96*se_ca, ca_hi=p_ca+1.96*se_ca,
       nc_est=p_nc, nc_lo=p_nc-1.96*se_nc, nc_hi=p_nc+1.96*se_nc,
       diff=diff, diff_lo=diff-1.96*se_diff, diff_hi=diff+1.96*se_diff)
}

calc_stratified <- function(model, design, stratum_var, stratum_levels) {
  design$variables$cc_fig <- complete.cases(design$variables[, vars_b])
  d_cc <- subset(design, cc_fig); df_cc <- d_cc$variables
  
  results <- lapply(stratum_levels, function(lv) {
    idx <- df_cc[[stratum_var]] == lv & !is.na(df_cc[[stratum_var]])
    n_ca <- sum(df_cc$cancer_history[idx]=="Cancer history",    na.rm=TRUE)
    n_nc <- sum(df_cc$cancer_history[idx]=="No cancer history", na.rm=TRUE)
    if (n_ca < 25 || n_nc < 25) {
      cat(sprintf("  %s: skipped (Cancer n=%d, No cancer n=%d)\n", lv, n_ca, n_nc))
      return(NULL)
    }
    pred_ca <- rep(NA_real_, nrow(df_cc))
    pred_nc <- rep(NA_real_, nrow(df_cc))
    pred_ca[idx] <- as.numeric(predict(model,
                                       newdata = df_cc[idx,] %>% mutate(cancer_history=factor("Cancer history",    levels=levels(analytic_df$cancer_history))),
                                       type="response"))
    pred_nc[idx] <- as.numeric(predict(model,
                                       newdata = df_cc[idx,] %>% mutate(cancer_history=factor("No cancer history", levels=levels(analytic_df$cancer_history))),
                                       type="response"))
    d_cc$variables$pred_ca_s <- pred_ca
    d_cc$variables$pred_nc_s <- pred_nc
    d_lv <- subset(d_cc, get(stratum_var) == lv)
    obj_ca <- svymean(~pred_ca_s, d_lv, na.rm=TRUE); obj_nc <- svymean(~pred_nc_s, d_lv, na.rm=TRUE)
    p_ca <- 100*as.numeric(coef(obj_ca)); se_ca <- 100*as.numeric(SE(obj_ca))
    p_nc <- 100*as.numeric(coef(obj_nc)); se_nc <- 100*as.numeric(SE(obj_nc))
    diff <- p_ca - p_nc; se_diff <- sqrt(se_ca^2 + se_nc^2)
    list(stratum=lv,
         ca_est=p_ca, ca_lo=p_ca-1.96*se_ca, ca_hi=p_ca+1.96*se_ca,
         nc_est=p_nc, nc_lo=p_nc-1.96*se_nc, nc_hi=p_nc+1.96*se_nc,
         diff=diff, diff_lo=diff-1.96*se_diff, diff_hi=diff+1.96*se_diff)
  })
  Filter(Negate(is.null), results)
}

results_to_df <- function(res_list, panel_label) {
  bind_rows(lapply(res_list, function(r) {
    data.frame(panel=panel_label, stratum=r$stratum,
               cancer_status=c("Cancer history","No cancer history"),
               est=c(r$ca_est, r$nc_est), lo=c(r$ca_lo, r$nc_lo), hi=c(r$ca_hi, r$nc_hi),
               diff_label=sprintf("%.1f pp (%.1f, %.1f)", r$diff, r$diff_lo, r$diff_hi),
               stringsAsFactors=FALSE)
  }))
}

# ==============================================================================
# 9-1  Panel A — Overall
# ==============================================================================
cat("--- Panel A: Overall ---\n")
res_a <- calc_overall(model_b, nhis_design)
cat(sprintf("  Cancer history:    %.1f%% (%.1f, %.1f)\n", res_a$ca_est, res_a$ca_lo, res_a$ca_hi))
cat(sprintf("  No cancer history: %.1f%% (%.1f, %.1f)\n", res_a$nc_est, res_a$nc_lo, res_a$nc_hi))
cat(sprintf("  Difference:        %.1f pp (%.1f, %.1f)\n\n", res_a$diff, res_a$diff_lo, res_a$diff_hi))

df_a <- data.frame(
  panel="A. Overall (Model B standardized)", stratum="Overall",
  cancer_status=c("Cancer history","No cancer history"),
  est=c(res_a$ca_est, res_a$nc_est), lo=c(res_a$ca_lo, res_a$nc_lo), hi=c(res_a$ca_hi, res_a$nc_hi),
  diff_label=sprintf("Difference: \u2212%.1f pp (95%% CI: %.1f, %.1f)",
                     abs(res_a$diff), res_a$diff_lo, res_a$diff_hi),
  stringsAsFactors=FALSE
)
plot_df <- df_a

# ==============================================================================
# 9-2  Panel B — By Insurance Type
# ==============================================================================
if (include_panel_b) {
  cat("--- Panel B: By Insurance Type ---\n")
  res_b_list <- calc_stratified(model_b, nhis_design, "insurance_type",
                                c("Private","Medicare","Medicaid"))
  if (length(res_b_list) > 0) {
    plot_df <- bind_rows(plot_df,
                         results_to_df(res_b_list, sprintf("B. By Insurance Type (P for interaction = %.3f)", p_int_ins)))
    cat("\n")
  }
}

# ==============================================================================
# 9-3  Panel C — By Race/Ethnicity
# ==============================================================================
if (include_panel_c) {
  cat("--- Panel C: By Race/Ethnicity ---\n")
  res_c_list <- calc_stratified(model_b, nhis_design, "race_ethnicity",
                                c("Non-Hispanic White","Non-Hispanic Black","Hispanic"))
  if (length(res_c_list) > 0) {
    plot_df <- bind_rows(plot_df,
                         results_to_df(res_c_list, sprintf("C. By Race/Ethnicity (P for interaction = %.3f)", p_int_race)))
    cat("\n")
  }
}

# ==============================================================================
# 9-4  Build Figure 1
# ==============================================================================
cat("--- Building Figure 1 ---\n")

plot_df$cancer_status <- factor(plot_df$cancer_status,
                                levels=c("Cancer history","No cancer history"))
plot_df$panel <- factor(plot_df$panel, levels=unique(plot_df$panel))
y_max <- ceiling(max(plot_df$hi, na.rm=TRUE) / 5) * 5 + 5

annot_df <- df_a %>%
  filter(cancer_status == "Cancer history") %>%
  select(panel, stratum, diff_label) %>% distinct()

fig1 <- ggplot(plot_df, aes(x=stratum, y=est, fill=cancer_status)) +
  geom_col(position=position_dodge(width=0.75), width=0.65) +
  geom_errorbar(aes(ymin=lo, ymax=hi),
                position=position_dodge(width=0.75), width=0.25, linewidth=0.5) +
  geom_text(data=annot_df, aes(x=stratum, y=y_max-1, label=diff_label),
            inherit.aes=FALSE, size=3.2, fontface="italic", hjust=0.5) +
  facet_wrap(~panel, scales="free_x", ncol=1) +
  scale_fill_manual(
    values = c("Cancer history"="#BDBDBD", "No cancer history"="#525252"),
    labels = c("Cancer history","No cancer history")
  ) +
  scale_y_continuous(limits=c(0, y_max), breaks=seq(0, y_max, by=5), expand=c(0,0)) +
  labs(
    title = "Figure 1. Adjusted Predicted Probability of Current GLP-1 RA Use\nby Cancer History (NHIS 2024)",
    x = NULL,
    y = "Adjusted Predicted Probability of GLP-1 RA Use (%)",
    fill = NULL,
    caption = paste(
      "Model B standardized estimates (marginal standardization). Error bars = 95% CI.",
      ifelse(include_panel_b, sprintf("Panel B: P for interaction (Cancer\u00d7Insurance) = %.3f.", p_int_ins), ""),
      ifelse(include_panel_c, sprintf("Panel C: P for interaction (Cancer\u00d7Race/Ethnicity) = %.3f.", p_int_race), "")
    )
  ) +
  theme_minimal(base_size=12) +
  theme(
    plot.title       = element_text(hjust=0.5, face="bold", size=13),
    plot.caption     = element_text(size=8, hjust=0, color="grey40"),
    legend.position  = "bottom",
    legend.key.size  = unit(0.5,"cm"),
    axis.text.x      = element_text(angle=30, hjust=1, size=10),
    axis.text.y      = element_text(size=10),
    axis.title.y     = element_text(size=11),
    strip.text       = element_text(face="bold", size=10),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

# ==============================================================================
# 9-5  Save
# ==============================================================================
print(fig1)
n_panels   <- length(unique(plot_df$panel))
fig_height <- 3.5 + 2.5 * (n_panels - 1)

ggsave("Figure1_adjusted_probabilities.png", fig1, width=9, height=fig_height, dpi=300)
ggsave("Figure1_adjusted_probabilities.pdf", fig1, width=9, height=fig_height)

cat("\n=============================================================================\n")
cat("Block 09 complete\n")
cat(sprintf("  Panel A: %.1f pp gap (%.1f, %.1f)\n", res_a$diff, res_a$diff_lo, res_a$diff_hi))
cat(sprintf("  Panel B: %s\n", ifelse(include_panel_b, "INCLUDED", "omitted (P >= 0.10)")))
cat(sprintf("  Panel C: %s\n", ifelse(include_panel_c, "INCLUDED", "omitted (P >= 0.10)")))
cat("  Saved: Figure1_adjusted_probabilities.png / .pdf\n")
cat("=============================================================================\n")


# ==============================================================================
# NHIS 2024 Analysis Pipeline — Block 10: Supplemental Race/Ethnicity Table
# ==============================================================================
# Input:  analytic_df (Block 01), nhis_design (Block 02),
#         model_b (Block 04), model_c (Block 08)
# Output: GLP1_supplemental_race_table.docx
#
# Produces two outputs:
#   1. eTable 3 key predictors (Model B vs Model C) — reproduced for reference
#   2. Race/ethnicity-stratified standardized estimates from Model B and C,
#      to assess whether attenuation varies by racial/ethnic subgroup
#
# Also tests P for interaction (cancer x race/ethnicity) in Model C.
# Cell suppression threshold: n >= 30 per NCHS guidelines.
# ==============================================================================

library(dplyr)
library(survey)
library(flextable)
library(officer)

if (!exists("nhis_design")) stop("nhis_design not found. Run Block 02 first.")
if (!exists("model_b"))     stop("model_b not found. Run Block 04 first.")
if (!exists("model_c"))     stop("model_c not found. Run Block 08 first.")

NCHS_THRESHOLD <- 30

cat("=============================================================================\n")
cat("Block 10: Supplemental Race/Ethnicity Table\n")
cat("=============================================================================\n\n")

# Sanity check
need_vars <- c("glp1_use","cancer_history","race_ethnicity","age_years","sex",
               "education","income_fpl_2cat","insurance_type","region",
               "bmi_category","ascvd","ckd","cost_barrier","spd_k6")
missing_v <- setdiff(need_vars, names(analytic_df))
if (length(missing_v) > 0) stop("Missing variables: ", paste(missing_v, collapse=", "))
cat("All required variables confirmed\n\n")

# ==============================================================================
# Formatting helpers
# ==============================================================================
extract_coef <- function(model, coef_name) {
  coefs <- coef(model); vcovs <- vcov(model); summ <- summary(model)$coefficients
  if (!coef_name %in% names(coefs)) return(list(or=NA, lo=NA, hi=NA, p=NA))
  est <- coefs[coef_name]; se <- sqrt(vcovs[coef_name, coef_name])
  list(or=exp(est), lo=exp(est-1.96*se), hi=exp(est+1.96*se), p=summ[coef_name,"Pr(>|t|)"])
}
fmt_or  <- function(res) { if (any(is.na(c(res$or,res$lo,res$hi)))) return("\u2014"); sprintf("%.2f (%.2f, %.2f)", res$or, res$lo, res$hi) }
fmt_p   <- function(p)   { if (is.null(p)||is.na(p)) return("\u2014"); if (p<0.001) return("<.001"); sprintf("%.3f", p) }
fmt_pct <- function(est, lo, hi) { if (any(is.na(c(est,lo,hi)))) return("\u2014"); sprintf("%.1f (%.1f, %.1f)", est, lo, hi) }

# ==============================================================================
# Section 1: eTable 3 key predictors (reuse models, no re-fitting)
# ==============================================================================
cat("--- Section 1: Model B vs Model C key predictors ---\n")
ca_b   <- extract_coef(model_b, "cancer_historyCancer history")
ca_c   <- extract_coef(model_c, "cancer_historyCancer history")
cost_c <- extract_coef(model_c, "cost_barrierYes")
spd_c  <- extract_coef(model_c, "spd_k6Yes")
pct_atten <- 100 * (log(ca_c$or) - log(ca_b$or)) / abs(log(ca_b$or))
cat(sprintf("  Model B cancer aOR: %.2f (%.2f, %.2f)\n", ca_b$or, ca_b$lo, ca_b$hi))
cat(sprintf("  Model C cancer aOR: %.2f (%.2f, %.2f)\n", ca_c$or, ca_c$lo, ca_c$hi))
cat(sprintf("  Attenuation: %.1f%%\n\n", pct_atten))

vars_b <- c("glp1_use","cancer_history","age_years","sex","race_ethnicity",
            "education","income_fpl_2cat","insurance_type","region","bmi_category","ascvd","ckd")
n_b <- sum(complete.cases(analytic_df[, vars_b]))
n_c <- sum(complete.cases(analytic_df[, c(vars_b,"cost_barrier","spd_k6")]))

et3_key <- data.frame(
  Predictor = c("Cancer history","Cost barrier","Serious psychological distress (K6)","Other covariates","Sample size"),
  Category  = c("Cancer vs No cancer history",
                "Any unmet or delayed care due to cost (Yes vs No)",
                "Kessler-6 score indicating serious psychological distress (Yes vs No)",
                "Age, sex, race/ethnicity, education, income, insurance, region, BMI, ASCVD, CKD",""),
  OR_B = c(fmt_or(ca_b),"\u2014","\u2014","Included (see eTable 1)",sprintf("n = %s", formatC(n_b,format="d",big.mark=","))),
  p_B  = c(fmt_p(ca_b$p),"\u2014","\u2014","\u2014","\u2014"),
  OR_C = c(fmt_or(ca_c),fmt_or(cost_c),fmt_or(spd_c),"Included (see eTable 1)",sprintf("n = %s", formatC(n_c,format="d",big.mark=","))),
  p_C  = c(fmt_p(ca_c$p),fmt_p(cost_c$p),fmt_p(spd_c$p),"\u2014","\u2014"),
  stringsAsFactors=FALSE
)
colnames(et3_key) <- c("Predictor","Category","Model B OR (95% CI)","Model B P value","Model C OR (95% CI)","Model C P value")

# ==============================================================================
# Section 2: Race interaction P-values
# ==============================================================================
cat("--- Section 2: Race interaction P-values ---\n")

safe_p_int <- function(model_int, term_pattern) {
  tryCatch({
    wt <- regTermTest(model_int, as.formula(paste0("~", term_pattern)), method="Wald")
    as.numeric(wt$p)
  }, error = function(e) {
    coefs <- summary(model_int)$coefficients
    rows  <- grep(term_pattern, rownames(coefs))
    if (length(rows) > 0) min(coefs[rows,"Pr(>|t|)"], na.rm=TRUE) else NA_real_
  })
}

m_b_race_int <- tryCatch(
  svyglm(glp1_use ~ cancer_history * race_ethnicity + age_years + sex +
           education + income_fpl_2cat + insurance_type + region + bmi_category + ascvd + ckd,
         design=nhis_design, family=quasibinomial()),
  error=function(e) { cat("  Model B race interaction error:", conditionMessage(e), "\n"); NULL }
)

m_c_race_int <- tryCatch(
  svyglm(glp1_use ~ cancer_history * race_ethnicity + age_years + sex +
           education + income_fpl_2cat + insurance_type + region + bmi_category + ascvd + ckd +
           cost_barrier + spd_k6,
         design=nhis_design, family=quasibinomial()),
  error=function(e) { cat("  Model C race interaction error:", conditionMessage(e), "\n"); NULL }
)

p_race_b <- if (!is.null(m_b_race_int)) safe_p_int(m_b_race_int, "cancer_history:race_ethnicity") else NA_real_
p_race_c <- if (!is.null(m_c_race_int)) safe_p_int(m_c_race_int, "cancer_history:race_ethnicity") else NA_real_
cat(sprintf("  P for interaction (Cancer x Race): Model B = %s | Model C = %s\n\n", fmt_p(p_race_b), fmt_p(p_race_c)))

# ==============================================================================
# Section 3: Race-stratified marginal standardization
# ==============================================================================
cat("--- Section 3: Race-stratified estimates ---\n\n")

std_race <- function(model, design, race_levels, model_label) {
  mv <- intersect(all.vars(formula(model)), names(design$variables))
  design$variables$cc_b10 <- complete.cases(design$variables[, mv])
  d_cc <- subset(design, cc_b10); df_cc <- d_cc$variables
  
  bind_rows(lapply(race_levels, function(rg) {
    idx_r <- df_cc$race_ethnicity == rg & !is.na(df_cc$race_ethnicity)
    n_ca  <- sum(df_cc$cancer_history[idx_r]=="Cancer history",    na.rm=TRUE)
    n_nc  <- sum(df_cc$cancer_history[idx_r]=="No cancer history", na.rm=TRUE)
    if (n_ca < NCHS_THRESHOLD || n_nc < NCHS_THRESHOLD) {
      cat(sprintf("  %s | %s: suppressed (Cancer n=%d, No cancer n=%d)\n", model_label, rg, n_ca, n_nc))
      return(data.frame(Model=model_label, Race=rg, Cancer_pct="\u2014", NoCancer_pct="\u2014", Gap="\u2014", stringsAsFactors=FALSE))
    }
    pred_ca <- rep(NA_real_, nrow(df_cc))
    pred_nc <- rep(NA_real_, nrow(df_cc))
    pred_ca[idx_r] <- as.numeric(predict(model,
                                         newdata=df_cc[idx_r,] %>% mutate(cancer_history=factor("Cancer history",    levels=levels(analytic_df$cancer_history))),
                                         type="response"))
    pred_nc[idx_r] <- as.numeric(predict(model,
                                         newdata=df_cc[idx_r,] %>% mutate(cancer_history=factor("No cancer history", levels=levels(analytic_df$cancer_history))),
                                         type="response"))
    d_cc$variables$pred_ca_r <- pred_ca
    d_cc$variables$pred_nc_r <- pred_nc
    d_r    <- subset(d_cc, race_ethnicity==rg)
    obj_ca <- svymean(~pred_ca_r, d_r, na.rm=TRUE); obj_nc <- svymean(~pred_nc_r, d_r, na.rm=TRUE)
    p_ca <- 100*as.numeric(coef(obj_ca)); se_ca <- 100*as.numeric(SE(obj_ca))
    p_nc <- 100*as.numeric(coef(obj_nc)); se_nc <- 100*as.numeric(SE(obj_nc))
    diff <- p_ca - p_nc; se_diff <- sqrt(se_ca^2 + se_nc^2)
    cat(sprintf("  %s | %s: Cancer %.1f%% | No cancer %.1f%% | Gap %.1f pp\n", model_label, rg, p_ca, p_nc, diff))
    data.frame(Model=model_label, Race=rg,
               Cancer_pct  = fmt_pct(p_ca, p_ca-1.96*se_ca, p_ca+1.96*se_ca),
               NoCancer_pct= fmt_pct(p_nc, p_nc-1.96*se_nc, p_nc+1.96*se_nc),
               Gap         = fmt_pct(diff, diff-1.96*se_diff, diff+1.96*se_diff),
               stringsAsFactors=FALSE)
  }))
}

race_levels_use <- c("Non-Hispanic White","Non-Hispanic Black","Hispanic")
race_b_df <- std_race(model_b, nhis_design, race_levels_use, "Model B"); cat("\n")
race_c_df <- std_race(model_c, nhis_design, race_levels_use, "Model C"); cat("\n")

race_table <- bind_rows(race_b_df, race_c_df)
colnames(race_table) <- c("Model","Race/Ethnicity",
                          "Cancer history: standardized % (95% CI)",
                          "No cancer history: standardized % (95% CI)",
                          "Absolute gap (Cancer \u2212 No cancer), pp (95% CI)")

# ==============================================================================
# Flextables
# ==============================================================================
ft_et3_key <- flextable(et3_key) %>%
  theme_booktabs() %>% autofit() %>%
  align(align="left",   j=1:2, part="all") %>%
  align(align="center", j=3:6, part="body") %>%
  align(align="center", j=3:6, part="header") %>%
  font(fontname="Times New Roman", part="all") %>%
  fontsize(size=10, part="body") %>% fontsize(size=10, part="header") %>%
  bold(part="header") %>% width(j=2, width=2.5)

ft_race <- flextable(race_table) %>%
  theme_booktabs() %>% autofit() %>%
  align(align="left",   j=1:2, part="all") %>%
  align(align="center", j=3:5, part="body") %>%
  align(align="center", j=3:5, part="header") %>%
  font(fontname="Times New Roman", part="all") %>%
  fontsize(size=10, part="body") %>% fontsize(size=10, part="header") %>%
  bold(part="header") %>% merge_v(j=1) %>% valign(j=1, valign="top")

footnote_et3_key <- paste(
  "Model B: adjusted for cancer history, age (continuous), sex, race/ethnicity (4 categories),",
  "education, income (<200% vs \u2265200% FPL), insurance (Private/Medicare/Medicaid),",
  "region, BMI (3 categories), ASCVD, and CKD.",
  "Model C: Model B + cost barrier + serious psychological distress (Kessler-6).",
  sprintf("Attenuation of cancer history OR (log-OR scale): %.1f%%.", pct_atten),
  "Model C is exploratory attenuation; cross-sectional design precludes causal mediation inference.",
  "OR = odds ratio; CI = confidence interval."
)

footnote_race <- paste(
  "Standardized prevalences computed via marginal standardization within each race/ethnicity stratum,",
  "survey-weighted. 'Other' race/ethnicity excluded due to heterogeneity and small cell sizes.",
  sprintf("Strata with fewer than %d unweighted observations per cancer history group are suppressed,", NCHS_THRESHOLD),
  "consistent with NCHS analytic guidelines for NHIS public-use data.",
  sprintf("P for interaction (Cancer \u00d7 Race/Ethnicity): Model B = %s; Model C = %s.", fmt_p(p_race_b), fmt_p(p_race_c)),
  "pp = percentage points; CI = confidence interval."
)

# Export
doc10 <- read_docx() %>%
  body_add_par("eTable 3. Exploratory Attenuation Model: Model B vs Model C Key Predictors (NHIS 2024)", style="heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft_et3_key) %>%
  body_add_par("") %>%
  body_add_par(footnote_et3_key, style="Normal") %>%
  body_add_par("") %>%
  body_add_par("Supplemental Table. Standardized GLP-1 RA Use by Cancer History and Race/Ethnicity: Model B vs Model C (NHIS 2024)", style="heading 1") %>%
  body_add_par("") %>%
  body_add_flextable(ft_race) %>%
  body_add_par("") %>%
  body_add_par(footnote_race, style="Normal")

print(doc10, target="GLP1_supplemental_race_table.docx")

cat("\n=============================================================================\n")
cat("Block 10 complete \u2014 GLP1_supplemental_race_table.docx saved\n")
cat(sprintf("  Cancer OR attenuation B->C: %.1f%%\n", pct_atten))
cat(sprintf("  P for interaction (Cancer x Race): Model B = %s | Model C = %s\n", fmt_p(p_race_b), fmt_p(p_race_c)))
cat("=============================================================================\n")

