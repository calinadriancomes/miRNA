# ===========================
# run_pipeline.R
# Da capo al fine pipeline
# ===========================

rm(list = ls())
options(stringsAsFactors = FALSE)

source("R/01_config.R")
source("R/02_io.R")
source("R/03_preprocess.R")
source("R/04_mirna_stats.R")
source("R/05_survival.R")
source("R/06_penalized_cox.R")
source("R/07_rsf.R")
source("R/08_reporting.R")

dir.create(PATH$out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PATH$out_dir, "figures"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PATH$out_dir, "tables"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(PATH$out_dir, "models"), showWarnings = FALSE, recursive = TRUE)

# 1) Load
raw <- load_all_data(PATH$smcp_xlsx, PATH$mirna_xlsx)

# 2) Harmonize / join
dat <- build_analysis_dataset(raw$smcp, raw$mirna, cfg = CFG)

# sanity check
message("Rows in analysis dataset: ", nrow(dat))
message("Columns: ", paste(colnames(dat), collapse = ", "))

# 3) miRNA expression analyses (1–18)
expr_results <- run_expression_suite(dat, cfg = CFG, out_dir = PATH$out_dir)

# 4) Survival analyses (19 + Cox)
surv_results <- run_survival_suite(dat, cfg = CFG, out_dir = PATH$out_dir)

# 5) Penalized Cox (LASSO)
lasso_results <- run_lasso_cox(dat, cfg = CFG, out_dir = PATH$out_dir)

# 6) Random survival forest
rsf_results <- run_rsf(dat, cfg = CFG, out_dir = PATH$out_dir)

# 7) Export master tables
export_master_outputs(expr_results, surv_results, lasso_results, rsf_results, out_dir = PATH$out_dir)

message("DONE. Outputs are in: ", PATH$out_dir)