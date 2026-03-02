# ===========================
# R/01_config.R
# Config + column mapping
# ===========================

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(openxlsx)
  library(survival)
  library(survminer)
  library(glmnet)
  library(randomForestSRC)
  library(broom)
  library(forcats)
})

PATH <- list(
  smcp_xlsx = "SMPC.xlsx",
  mirna_xlsx = "miRNA_0.xlsx",
  out_dir   = "OUTPUT"
)

# ---- Column mapping ----
# Here you put the real names from Excel if they differ.
CFG <- list(
  # key columns
  id_col          = "ID",          # if you don't have an ID, use another common identifier
  diagnosis_col   = "Diagnosis",   # PV / ET / PMF / Control
  sex_col         = "Gender",      # Male/Female sau M/F
  age_col         = "AgeAtDiagnosis",
  # miRNA
  mir_125_col     = "miRNA-125b-5p_MeanDCq",
  mir_223_col     = "miRNA-223-3p_MeanDCq",

  # survival
  time_col        = "OS_months",   # follow-up time in months
  event_col       = "Event",       # 1=death, 0=censored

  # mutations
  jak2_col        = "JAK2",        # Positive/Negative or 1/0
  calr_col        = "CALR",
  mpl_col         = "MPL",

  # clinical binary covariates / columns mentioned by you (letters in Excel)
  major_thrombosis_col   = "MajorTrombosis",          # N
  arterial_thrombosis_col= "ArterialThrombosis",      # M
  venous_thrombosis_col  = "VenousThrombosis",        # O
  palpable_spleno_col    = "PalpableSplenomegaly",    # Q
  secondary_mf_col       = "SecondaryMielofibrosis",  # S
  risk_factors_col       = "RiskFactors",             # AD
  hydroxyurea_col        = "HydroxyureaTreatment",    # AL
  antiplatelet_col       = "TreatmentWithAntiplateletAgents", # AM
  anagrelide_col         = "TratamentWithAnagrelide", # AN
  interferon_col         = "InterferonTreatment",     # AQ
  blood_emissions_col    = "BloodEmissions",          # AO (phlebotomy)
  smoking_col            = "SmokingHabits",           # AS
  alcohol_col            = "AlcoholHabits",           # AU

  # lab thresholds for univ/multiv Cox set 22
  hgb_col   = "Hemoglobin",
  hct_col   = "Hematocrit",
  plt_col   = "Platelets",
  wbc_col   = "WBC",
  ldh_col   = "LDH"
)

# Diagnosis order for plotting
DIAG_LEVELS <- c("Control", "PV", "ET", "PMF")