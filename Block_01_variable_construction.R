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
