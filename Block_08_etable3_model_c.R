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
