# ===========================
# R/08_reporting.R
# Gather + export master files
# ===========================

export_master_outputs <- function(expr, surv, lasso, rsf, out_dir) {
  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "Log-rank")
  openxlsx::writeData(wb, "Log-rank", data.frame(
    chisq = surv$logrank$stat,
    df = surv$logrank$df,
    p = surv$logrank$p
  ))

  openxlsx::addWorksheet(wb, "Cox_opt")
  openxlsx::writeData(wb, "Cox_opt", surv$cox_opt_table)

  openxlsx::addWorksheet(wb, "Cox_univariable")
  openxlsx::writeData(wb, "Cox_univariable", surv$cox_uni_table)

  openxlsx::addWorksheet(wb, "Cox_multivariable")
  openxlsx::writeData(wb, "Cox_multivariable", surv$cox_multi_table)

  openxlsx::addWorksheet(wb, "PH_test_opt")
  openxlsx::writeData(wb, "PH_test_opt", surv$ph_test$table)

  openxlsx::addWorksheet(wb, "LASSO_selected_1se")
  openxlsx::writeData(wb, "LASSO_selected_1se", data.frame(variable = lasso$selected_1se))

  openxlsx::addWorksheet(wb, "LASSO_selected_min")
  openxlsx::writeData(wb, "LASSO_selected_min", data.frame(variable = lasso$selected_min))

  openxlsx::addWorksheet(wb, "RSF_importance")
  openxlsx::writeData(wb, "RSF_importance", rsf$importance)

  openxlsx::saveWorkbook(wb, file.path(out_dir, "tables", "MASTER_results.xlsx"), overwrite = TRUE)
}