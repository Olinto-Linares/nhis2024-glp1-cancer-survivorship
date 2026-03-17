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
