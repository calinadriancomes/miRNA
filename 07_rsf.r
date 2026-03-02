# ===========================
# R/07_rsf.R
# Random Survival Forest with factor conversion
# ===========================

run_rsf <- function(df, cfg, out_dir) {
  diag <- cfg$diagnosis_col
  time <- cfg$time_col
  event <- cfg$event_col

  pat <- df %>% filter(.data[[diag]] %in% c("PV","ET","PMF")) %>%
    filter(!is.na(.data[[time]]) & !is.na(.data[[event]]))

  pat$Male <- as.integer(pat[[cfg$sex_col]] == "Male")

  # Select RSF predictors (as in your example)
  rsf_vars <- c(
    cfg$diagnosis_col, "Age_gt60", "Male", "Hgb_lt16", "Hct_lt48", "Plt_lt450",
    "WBC_gt11", "JAK2_pos", "CALR_pos", "MPL_pos", cfg$ldh_col
  )

  # Convert characters to factors (fix the exact error you got)
  for (v in rsf_vars) {
    if (v %in% colnames(pat) && is.character(pat[[v]])) pat[[v]] <- factor(pat[[v]])
  }
  # ensure diagnosis factor
  if (cfg$diagnosis_col %in% colnames(pat)) pat[[cfg$diagnosis_col]] <- factor(pat[[cfg$diagnosis_col]])

  # Build formula
  f <- as.formula(
    paste0("Surv(", time, ",", event, ") ~ ", paste(rsf_vars, collapse = " + "))
  )

  set.seed(1)
  rsf <- rfsrc(f, data = pat, ntree = 2000, importance = TRUE)

  # Importance
  imp <- rsf$importance
  imp_df <- data.frame(variable = names(imp), importance = as.numeric(imp), row.names = NULL) %>%
    arrange(desc(importance))

  openxlsx::write.xlsx(imp_df, file.path(out_dir, "tables", "rsf_variable_importance.xlsx"),
                       overwrite = TRUE)

  # quick importance plot
  p <- ggplot(imp_df, aes(x = reorder(variable, importance), y = importance)) +
    geom_col() +
    coord_flip() +
    labs(x = "", y = "Variable importance", title = "Random Survival Forest: variable importance") +
    theme_minimal(base_size = 12)
  ggsave(file.path(out_dir, "figures", "rsf_variable_importance.png"), p, width = 7.2, height = 5.2, dpi = 300)

  list(rsf = rsf, importance = imp_df, importance_plot = p)
}