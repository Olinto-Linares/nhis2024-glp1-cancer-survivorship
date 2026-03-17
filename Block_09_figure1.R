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
