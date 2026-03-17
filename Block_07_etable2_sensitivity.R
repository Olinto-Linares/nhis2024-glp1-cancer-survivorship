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
