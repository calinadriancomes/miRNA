# ===========================
# R/05_survival.R
# KM + log-rank + Cox + PH tests + forest plot + tables
# ===========================

km_plot <- function(df, cfg, title, out_png) {
  time <- cfg$time_col
  event <- cfg$event_col
  diag <- cfg$diagnosis_col

  fit <- survfit(Surv(.data[[time]], .data[[event]]) ~ .data[[diag]], data = df)

  p <- ggsurvplot(
    fit, data = df, risk.table = TRUE, pval = FALSE,
    conf.int = FALSE, legend.title = "Diagnosis",
    xlab = "Months", ylab = "Overall survival probability",
    title = title
  )

  ggsave(out_png, p$plot, width = 7.8, height = 5.2, dpi = 300)
  # risk table optional export: only plot saved; can save p$table similarly
  list(fit = fit, plot = p)
}

logrank_test <- function(df, cfg) {
  time <- cfg$time_col
  event <- cfg$event_col
  diag <- cfg$diagnosis_col

  lr <- survdiff(Surv(df[[time]], df[[event]]) ~ df[[diag]])
  p <- 1 - pchisq(lr$chisq, df = (length(lr$n) - 1))
  list(stat = lr$chisq, df = length(lr$n) - 1, p = p)
}

cox_pub_table <- function(fit) {
  s <- summary(fit)
  co <- as.data.frame(s$coefficients)
  ci <- as.data.frame(s$conf.int)

  out <- data.frame(
    term = rownames(co),
    HR = ci$`exp(coef)`,
    CI_low = ci$`lower .95`,
    CI_high = ci$`upper .95`,
    p = co$`Pr(>|z|)`,
    stringsAsFactors = FALSE
  )
  out
}

forest_plot_cox <- function(fit, title, out_png) {
  tb <- cox_pub_table(fit)
  tb <- tb %>%
    mutate(term = fct_rev(factor(term))) %>%
    mutate(logHR = log(HR),
           lo = log(CI_low),
           hi = log(CI_high))

  p <- ggplot(tb, aes(x = logHR, y = term)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point() +
    geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0.2) +
    labs(x = "log(HR)", y = "", title = title) +
    theme_minimal(base_size = 12)

  ggsave(out_png, p, width = 7.4, height = 4.6, dpi = 300)
  p
}

test_ph_assumption <- function(fit) {
  z <- cox.zph(fit)
  list(zph = z, table = as.data.frame(z$table))
}

run_survival_suite <- function(df, cfg, out_dir) {
  diag <- cfg$diagnosis_col
  time <- cfg$time_col
  event <- cfg$event_col

  # restrict to PV/ET/PMF (patients)
  pat <- df %>% filter(.data[[diag]] %in% c("PV","ET","PMF"))
  pat <- pat %>% filter(!is.na(.data[[time]]) & !is.na(.data[[event]]) & !is.na(.data[[diag]]))

  # KM separate per diagnosis
  km_sep <- list()
  for (d in c("PV","ET","PMF")) {
    sub <- pat %>% filter(.data[[diag]] == d)
    km_sep[[d]] <- km_plot(sub %>% mutate(One=d), # dummy group
                           cfg = modifyList(cfg, list(diagnosis_col="One")),
                           title = paste0("Overall survival: ", d),
                           out_png = file.path(out_dir,"figures",paste0("KM_",d,".png")))
  }

  # Combined KM PV-ET-PMF
  km_comb <- km_plot(pat, cfg, "Overall survival by diagnosis (PV vs ET vs PMF)",
                     file.path(out_dir,"figures","KM_PV_ET_PMF_combined.png"))

  # log-rank among the 3 groups
  lr <- logrank_test(pat, cfg)

  # Cox: Diagnosis + Age + JAK2 (optimized model you requested)
  pat$Diagnosis <- pat[[diag]]
  pat$Age <- pat[[cfg$age_col]]
  pat$JAK2_pos <- pat$JAK2_pos

  cox_opt <- coxph(Surv(pat[[time]], pat[[event]]) ~ Diagnosis + Age + JAK2_pos, data = pat)
  cox_opt_tab <- cox_pub_table(cox_opt)
  openxlsx::write.xlsx(cox_opt_tab,
                       file.path(out_dir,"tables","cox_publication_table_Diagnosis_Age_JAK2.xlsx"),
                       overwrite = TRUE)
  forest_plot_cox(cox_opt, "Cox model: Diagnosis + Age + JAK2",
                  file.path(out_dir,"figures","cox_forest_Diagnosis_Age_JAK2.png"))

  # PH test
  ph <- test_ph_assumption(cox_opt)
  openxlsx::write.xlsx(ph$table,
                       file.path(out_dir,"tables","cox_PH_test_Diagnosis_Age_JAK2.xlsx"),
                       overwrite = TRUE)

  # (22) Univariable Cox models for listed variables
  # Variables per your list:
  # Age_gt60, Male, Hgb_lt16, Hct_lt48, Plt_lt450, WBC_gt11, JAK2_pos, CALR_pos, LDH (+ MPL if desired)
  pat$Male <- as.integer(pat[[cfg$sex_col]] == "Male")
  pat$CALR_pos <- pat$CALR_pos
  pat$MPL_pos <- pat$MPL_pos
  pat$LDH <- pat[[cfg$ldh_col]]

  uni_vars <- c("Age_gt60","Male","Hgb_lt16","Hct_lt48","Plt_lt450","WBC_gt11","JAK2_pos","CALR_pos","MPL_pos","LDH")
  uni_tables <- lapply(uni_vars, function(v) {
    f <- as.formula(paste0("Surv(", time, ",", event, ") ~ ", v))
    fit <- coxph(f, data = pat)
    tb <- cox_pub_table(fit)
    tb$variable <- v
    tb
  }) %>% bind_rows()

  openxlsx::write.xlsx(uni_tables, file.path(out_dir,"tables","cox_univariable_publication_table.xlsx"),
                       overwrite = TRUE)

  # Multivariable Cox model (include the set)
  # You can adjust formula here if you want to keep Diagnosis always in.
  multiv_formula <- as.formula(
    paste0("Surv(", time, ",", event, ") ~ Diagnosis + Age_gt60 + Male + Hgb_lt16 + Hct_lt48 + Plt_lt450 + WBC_gt11 + JAK2_pos + CALR_pos + MPL_pos + LDH")
  )
  cox_multi <- coxph(multiv_formula, data = pat)
  cox_multi_tab <- cox_pub_table(cox_multi)
  openxlsx::write.xlsx(cox_multi_tab, file.path(out_dir,"tables","cox_multivariable_publication_table.xlsx"),
                       overwrite = TRUE)
  forest_plot_cox(cox_multi, "Multivariable Cox model",
                  file.path(out_dir,"figures","cox_multivariable_forest.png"))

  list(
    km_sep = km_sep,
    km_combined = km_comb,
    logrank = lr,
    cox_opt = cox_opt,
    cox_opt_table = cox_opt_tab,
    ph_test = ph,
    cox_uni_table = uni_tables,
    cox_multi = cox_multi,
    cox_multi_table = cox_multi_tab
  )
}