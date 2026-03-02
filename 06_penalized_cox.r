# ===========================
# R/06_penalized_cox.R
# LASSO Cox (glmnet) with fix for non-positive times
# ===========================

prep_glmnet_matrix <- function(df, vars, drop_na = TRUE) {
  d <- df[, vars, drop = FALSE]
  if (drop_na) d <- d[complete.cases(d), , drop = FALSE]

  # convert characters to factors, factors kept
  for (j in seq_along(d)) {
    if (is.character(d[[j]])) d[[j]] <- factor(d[[j]])
  }

  mm <- model.matrix(~ . , data = d)[, -1, drop = FALSE] # remove intercept
  list(mm = mm, df = d)
}

run_lasso_cox <- function(df, cfg, out_dir) {
  diag <- cfg$diagnosis_col
  time <- cfg$time_col
  event <- cfg$event_col

  pat <- df %>% filter(.data[[diag]] %in% c("PV","ET","PMF")) %>%
    filter(!is.na(.data[[time]]) & !is.na(.data[[event]]))

  # Ensure strictly positive times already by fix_survival_time()
  pat[[time]] <- fix_survival_time(pat[[time]])

  # Candidate predictors (you can expand)
  cand <- c(cfg$diagnosis_col, cfg$age_col, "Age_gt60",
            cfg$sex_col, "Male",
            "Hgb_lt16","Hct_lt48","Plt_lt450","WBC_gt11",
            "JAK2_pos","CALR_pos","MPL_pos", cfg$ldh_col)

  # Create Male for modeling
  pat$Male <- as.integer(pat[[cfg$sex_col]] == "Male")

  # Build X, y
  dat_mod <- pat %>% select(all_of(cand), all_of(time), all_of(event))

  # Convert diagnosis/sex to factors
  if (is.character(dat_mod[[cfg$diagnosis_col]])) dat_mod[[cfg$diagnosis_col]] <- factor(dat_mod[[cfg$diagnosis_col]])
  if (is.character(dat_mod[[cfg$sex_col]])) dat_mod[[cfg$sex_col]] <- factor(dat_mod[[cfg$sex_col]])

  # Replace Age col numeric
  dat_mod[[cfg$age_col]] <- as.numeric(dat_mod[[cfg$age_col]])
  dat_mod[[cfg$ldh_col]] <- as.numeric(dat_mod[[cfg$ldh_col]])

  # Prepare matrix
  Xprep <- prep_glmnet_matrix(dat_mod %>% select(-all_of(c(time, event))), vars = setdiff(names(dat_mod), c(time,event)), drop_na = TRUE)
  X <- Xprep$mm
  keep_idx <- rownames(X) # rows kept after complete.cases

  y <- Surv(dat_mod[as.integer(keep_idx), time], dat_mod[as.integer(keep_idx), event])

  set.seed(1)
  cvfit <- cv.glmnet(X, y, family = "cox", alpha = 1)

  # Selected coefficients at lambda.1se and lambda.min
  beta_1se <- coef(cvfit, s = "lambda.1se")
  beta_min <- coef(cvfit, s = "lambda.min")

  sel_1se <- rownames(beta_1se)[as.vector(beta_1se != 0)]
  sel_min <- rownames(beta_min)[as.vector(beta_min != 0)]

  out <- list(
    cvfit = cvfit,
    selected_1se = sel_1se,
    selected_min = sel_min
  )

  # export selections
  openxlsx::write.xlsx(
    list(
      selected_lambda_1se = data.frame(variable = sel_1se),
      selected_lambda_min = data.frame(variable = sel_min)
    ),
    file.path(out_dir, "tables", "lasso_cox_selected_variables.xlsx"),
    overwrite = TRUE
  )

  out
}