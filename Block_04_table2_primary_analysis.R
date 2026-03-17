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
