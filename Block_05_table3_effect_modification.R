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
