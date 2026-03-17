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
