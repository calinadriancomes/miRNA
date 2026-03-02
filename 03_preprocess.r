# ===========================
# R/03_preprocess.R
# Harmonize, recode, derived vars
# ===========================

as_binary01 <- function(x) {
  # robust conversion to 0/1 with common encodings
  if (is.logical(x)) return(as.integer(x))
  if (is.numeric(x)) return(as.integer(x > 0))
  z <- tolower(trimws(as.character(x)))
  as.integer(z %in% c("1","yes","y","true","positive","pos","mut","mutated"))
}

as_factor_posneg <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  z <- tolower(trimws(as.character(x)))
  out <- ifelse(z %in% c("1","yes","y","true","positive","pos","mut","mutated"), "Positive",
                ifelse(z %in% c("0","no","n","false","negative","neg","wt","wildtype"), "Negative", NA))
  factor(out, levels = c("Negative","Positive"))
}

normalize_sex <- function(x) {
  z <- tolower(trimws(as.character(x)))
  out <- ifelse(z %in% c("m","male","man","barbat","bărbat"), "Male",
                ifelse(z %in% c("f","female","woman","femeie"), "Female", NA))
  factor(out, levels = c("Female","Male"))
}

make_age_groups <- function(age) {
  cut(
    age,
    breaks = c(-Inf, 29, 49, 69, Inf),
    labels = c("18-29","30-49","50-69",">=70"),
    right = TRUE
  )
}

fix_survival_time <- function(t) {
  # glmnet Cox cannot take non-positive times; set 0 to a small epsilon
  t <- as.numeric(t)
  eps <- 1e-3
  t[is.na(t)] <- NA_real_
  t[t <= 0] <- eps
  t
}

build_analysis_dataset <- function(smcp_df, mirna_df, cfg) {
  stopifnot(cfg$id_col %in% colnames(smcp_df) || cfg$id_col %in% colnames(mirna_df))

  # join: prefer common ID; if missing, user must set cfg$id_col to a shared key
  df <- smcp_df %>%
    left_join(mirna_df, by = cfg$id_col)

  # diagnosis factor
  df[[cfg$diagnosis_col]] <- factor(as.character(df[[cfg$diagnosis_col]]), levels = DIAG_LEVELS)

  # sex
  df[[cfg$sex_col]] <- normalize_sex(df[[cfg$sex_col]])

  # age
  df[[cfg$age_col]] <- as.numeric(df[[cfg$age_col]])
  df$Age_group <- make_age_groups(df[[cfg$age_col]])
  df$Age_gt60  <- as.integer(df[[cfg$age_col]] > 60)

  # miRNA: ensure numeric
  df[[cfg$mir_125_col]] <- as.numeric(df[[cfg$mir_125_col]])
  df[[cfg$mir_223_col]] <- as.numeric(df[[cfg$mir_223_col]])

  # mutations
  df$JAK2_pos <- as_binary01(df[[cfg$jak2_col]])
  df$CALR_pos <- as_binary01(df[[cfg$calr_col]])
  df$MPL_pos  <- as_binary01(df[[cfg$mpl_col]])

  df$AnyMutation <- as.integer((df$JAK2_pos + df$CALR_pos + df$MPL_pos) > 0)
  df$TripleNegative <- as.integer((df$JAK2_pos + df$CALR_pos + df$MPL_pos) == 0)

  # survival
  df[[cfg$time_col]]  <- fix_survival_time(df[[cfg$time_col]])
  df[[cfg$event_col]] <- as_binary01(df[[cfg$event_col]])

  # labs and thresholds for set 22
  df[[cfg$hgb_col]] <- as.numeric(df[[cfg$hgb_col]])
  df[[cfg$hct_col]] <- as.numeric(df[[cfg$hct_col]])
  df[[cfg$plt_col]] <- as.numeric(df[[cfg$plt_col]])
  df[[cfg$wbc_col]] <- as.numeric(df[[cfg$wbc_col]])
  df[[cfg$ldh_col]] <- as.numeric(df[[cfg$ldh_col]])

  df$Hgb_lt16  <- as.integer(df[[cfg$hgb_col]] < 16)
  df$Hct_lt48  <- as.integer(df[[cfg$hct_col]] < 48)
  df$Plt_lt450 <- as.integer(df[[cfg$plt_col]] < 450000)
  df$WBC_gt11  <- as.integer(df[[cfg$wbc_col]] > 11e9) # adjust if units differ

  # clinical binary fields (coerce robustly)
  bin_cols <- c(
    cfg$major_thrombosis_col, cfg$arterial_thrombosis_col, cfg$venous_thrombosis_col,
    cfg$palpable_spleno_col, cfg$secondary_mf_col, cfg$risk_factors_col,
    cfg$hydroxyurea_col, cfg$antiplatelet_col, cfg$anagrelide_col,
    cfg$interferon_col, cfg$blood_emissions_col, cfg$smoking_col, cfg$alcohol_col
  )
  for (cc in bin_cols) {
    if (cc %in% colnames(df)) df[[cc]] <- as_binary01(df[[cc]])
  }

  df
}