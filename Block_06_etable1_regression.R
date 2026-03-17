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
