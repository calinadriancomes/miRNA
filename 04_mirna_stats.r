# ===========================
# R/04_mirna_stats.R
# Expression analyses (1–18)
# ===========================

wilcox_or_ttest <- function(x, g) {
  # x numeric, g factor with 2 levels
  x <- x[!is.na(x) & !is.na(g)]
  g <- droplevels(g[!is.na(x) & !is.na(g)])
  if (nlevels(g) != 2) return(list(test = NA, p = NA_real_, method = NA))
  # normality quick: if both groups >= 10 do Shapiro per group; else default Wilcoxon
  nx <- table(g)
  use_t <- FALSE
  if (all(nx >= 10)) {
    p1 <- tryCatch(shapiro.test(x[g == levels(g)[1]])$p.value, error = function(e) NA_real_)
    p2 <- tryCatch(shapiro.test(x[g == levels(g)[2]])$p.value, error = function(e) NA_real_)
    use_t <- !is.na(p1) && !is.na(p2) && (p1 > 0.05) && (p2 > 0.05)
  }
  if (use_t) {
    tt <- t.test(x ~ g)
    list(test = "t-test", p = tt$p.value, method = tt$method)
  } else {
    wt <- wilcox.test(x ~ g, exact = FALSE)
    list(test = "Mann–Whitney", p = wt$p.value, method = wt$method)
  }
}

summary_group <- function(df, xcol, bycols) {
  df %>%
    group_by(across(all_of(bycols))) %>%
    summarise(
      n = sum(!is.na(.data[[xcol]])),
      median = median(.data[[xcol]], na.rm = TRUE),
      q1 = quantile(.data[[xcol]], 0.25, na.rm = TRUE),
      q3 = quantile(.data[[xcol]], 0.75, na.rm = TRUE),
      mean = mean(.data[[xcol]], na.rm = TRUE),
      sd = sd(.data[[xcol]], na.rm = TRUE),
      .groups = "drop"
    )
}

make_boxplot <- function(df, x, group, title, out_png) {
  p <- ggplot(df, aes(x = .data[[group]], y = .data[[x]])) +
    geom_boxplot(outlier.alpha = 0.35) +
    geom_jitter(width = 0.12, height = 0, alpha = 0.35) +
    labs(x = group, y = x, title = title) +
    theme_minimal(base_size = 12)
  ggsave(out_png, p, width = 6.5, height = 4.2, dpi = 300)
  p
}

run_expression_suite <- function(df, cfg, out_dir) {
  diag <- cfg$diagnosis_col
  sex  <- cfg$sex_col
  ageg <- "Age_group"

  mirs <- c(cfg$mir_125_col, cfg$mir_223_col)

  tables <- list()
  plots  <- list()

  # Helper subsets
  patients <- df %>% filter(.data[[diag]] %in% c("PV","ET","PMF"))
  controls <- df %>% filter(.data[[diag]] == "Control")

  # (1) Patients vs Controls
  for (m in mirs) {
    sub <- df %>% filter(.data[[diag]] %in% c("Control","PV","ET","PMF")) %>%
      mutate(Group_PC = ifelse(.data[[diag]] == "Control", "Control", "Patient") %>% factor(levels=c("Control","Patient")))
    res <- wilcox_or_ttest(sub[[m]], sub$Group_PC)
    tables[[paste0("01_Patients_vs_Controls_", m)]] <- list(
      summary = summary_group(sub, m, "Group_PC"),
      test = data.frame(contrast="Patients vs Controls", miRNA=m, test=res$test, p=res$p)
    )
    plots[[paste0("box_01_", m)]] <- make_boxplot(sub, m, "Group_PC",
      title = paste0(m, " (MeanDCq): Patients vs Controls"),
      out_png = file.path(out_dir, "figures", paste0("box_01_Patients_vs_Controls_", m, ".png"))
    )
  }

  # (2) Patients vs Controls stratified by sex
  for (m in mirs) {
    sub <- df %>%
      filter(.data[[diag]] %in% c("Control","PV","ET","PMF"), !is.na(.data[[sex]])) %>%
      mutate(Group_PC = ifelse(.data[[diag]] == "Control", "Control", "Patient") %>% factor(levels=c("Control","Patient")))

    tab <- sub %>%
      group_by(.data[[sex]]) %>%
      group_modify(~{
        res <- wilcox_or_ttest(.x[[m]], .x$Group_PC)
        data.frame(sex = unique(.x[[sex]]), miRNA=m, test=res$test, p=res$p)
      }) %>% ungroup()

    tables[[paste0("02_PC_by_sex_", m)]] <- list(
      summary = summary_group(sub, m, c(sex, "Group_PC")),
      test = tab
    )

    # faceted boxplot
    p <- ggplot(sub, aes(x = Group_PC, y = .data[[m]])) +
      geom_boxplot(outlier.alpha = 0.35) +
      geom_jitter(width = 0.12, alpha = 0.35) +
      facet_wrap(vars(.data[[sex]])) +
      labs(x = "", y = m, title = paste0(m, " (MeanDCq): Patients vs Controls by sex")) +
      theme_minimal(base_size = 12)
    ggsave(file.path(out_dir,"figures",paste0("box_02_PC_by_sex_",m,".png")), p, width=8, height=4.2, dpi=300)
    plots[[paste0("box_02_", m)]] <- p
  }

  # (3) Within patients: PV vs ET vs PMF (Kruskal–Wallis)
  for (m in mirs) {
    sub <- patients %>% filter(!is.na(.data[[diag]]))
    kw <- kruskal.test(sub[[m]] ~ sub[[diag]])
    tables[[paste0("03_PV_ET_PMF_", m)]] <- list(
      summary = summary_group(sub, m, diag),
      test = data.frame(test="Kruskal–Wallis", miRNA=m, p=kw$p.value)
    )
    plots[[paste0("box_03_", m)]] <- make_boxplot(sub, m, diag,
      title = paste0(m, " (MeanDCq): PV vs ET vs PMF"),
      out_png = file.path(out_dir, "figures", paste0("box_03_PV_ET_PMF_", m, ".png"))
    )
  }

  # (4) PV/ET/PMF by sex (within each diagnosis: Male vs Female)
  for (m in mirs) {
    sub <- patients %>% filter(!is.na(.data[[sex]]))
    tab <- sub %>%
      group_by(.data[[diag]]) %>%
      group_modify(~{
        if (n_distinct(.x[[sex]], na.rm=TRUE) < 2) {
          data.frame(Diagnosis=unique(.x[[diag]]), miRNA=m, test=NA, p=NA_real_)
        } else {
          res <- wilcox_or_ttest(.x[[m]], .x[[sex]])
          data.frame(Diagnosis=unique(.x[[diag]]), miRNA=m, test=res$test, p=res$p)
        }
      }) %>% ungroup()
    tables[[paste0("04_Dx_by_sex_", m)]] <- list(
      summary = summary_group(sub, m, c(diag, sex)),
      test = tab
    )
  }

  # (5) PV/ET/PMF by age groups + Age>60
  for (m in mirs) {
    sub <- patients %>% filter(!is.na(.data[[ageg]]))
    # Age groups (KW per diagnosis)
    tab_ageg <- sub %>%
      group_by(.data[[diag]]) %>%
      group_modify(~{
        if (n_distinct(.x[[ageg]], na.rm=TRUE) < 2) {
          data.frame(Diagnosis=unique(.x[[diag]]), miRNA=m, test=NA, p=NA_real_)
        } else {
          kt <- kruskal.test(.x[[m]] ~ .x[[ageg]])
          data.frame(Diagnosis=unique(.x[[diag]]), miRNA=m, test="Kruskal–Wallis (Age groups)", p=kt$p.value)
        }
      }) %>% ungroup()

    # Age>60 (MW within diagnosis)
    sub$Age_gt60_fac <- factor(ifelse(sub$Age_gt60==1, ">60", "<=60"), levels=c("<=60",">60"))
    tab_gt60 <- sub %>%
      group_by(.data[[diag]]) %>%
      group_modify(~{
        if (n_distinct(.x$Age_gt60_fac, na.rm=TRUE) < 2) {
          data.frame(Diagnosis=unique(.x[[diag]]), miRNA=m, test=NA, p=NA_real_)
        } else {
          res <- wilcox_or_ttest(.x[[m]], .x$Age_gt60_fac)
          data.frame(Diagnosis=unique(.x[[diag]]), miRNA=m, test=res$test, p=res$p)
        }
      }) %>% ungroup()

    tables[[paste0("05_Dx_by_age_", m)]] <- list(
      summary_agegroups = summary_group(sub, m, c(diag, ageg)),
      test_agegroups = tab_ageg,
      summary_gt60 = summary_group(sub, m, c(diag, "Age_gt60_fac")),
      test_gt60 = tab_gt60
    )
  }

  # (6) Any mutation vs TripleNegative (within patients)
  for (m in mirs) {
    sub <- patients %>% mutate(MutGroup = factor(ifelse(AnyMutation==1, "AnyMutation", "TripleNegative"),
                                                levels=c("TripleNegative","AnyMutation")))
    res <- wilcox_or_ttest(sub[[m]], sub$MutGroup)
    tables[[paste0("06_Mutation_vs_TripleNeg_", m)]] <- list(
      summary = summary_group(sub, m, "MutGroup"),
      test = data.frame(contrast="Any mutation vs TripleNegative", miRNA=m, test=res$test, p=res$p)
    )
    plots[[paste0("box_06_", m)]] <- make_boxplot(sub, m, "MutGroup",
      title = paste0(m, " (MeanDCq): Any mutation vs Triple-negative"),
      out_png = file.path(out_dir, "figures", paste0("box_06_Mutation_vs_TripleNeg_", m, ".png"))
    )
  }

  # (7–18) Clinical/treatment factors within PV/ET/PMF
  # generic runner: for each factor run within each diagnosis comparing 0 vs 1
  factors_map <- list(
    MajorTrombosis = cfg$major_thrombosis_col,
    ArterialThrombosis = cfg$arterial_thrombosis_col,
    VenousThrombosis = cfg$venous_thrombosis_col,
    PalpableSplenomegaly = cfg$palpable_spleno_col,
    SecondaryMielofibrosis = cfg$secondary_mf_col,
    RiskFactors = cfg$risk_factors_col,
    HydroxyureaTreatment = cfg$hydroxyurea_col,
    TreatmentWithAntiplateletAgents = cfg$antiplatelet_col,
    TratamentWithAnagrelide = cfg$anagrelide_col,
    InterferonTreatment = cfg$interferon_col,
    BloodEmissions = cfg$blood_emissions_col,
    SmokingHabits = cfg$smoking_col,
    AlcoholHabits = cfg$alcohol_col
  )

  for (nm in names(factors_map)) {
    fcol <- factors_map[[nm]]
    if (!fcol %in% colnames(patients)) next

    for (m in mirs) {
      sub <- patients %>% filter(!is.na(.data[[fcol]]))
      sub$F <- factor(ifelse(sub[[fcol]]==1, "Yes", "No"), levels=c("No","Yes"))

      tab <- sub %>%
        group_by(.data[[diag]]) %>%
        group_modify(~{
          if (n_distinct(.x$F, na.rm=TRUE) < 2) {
            data.frame(Diagnosis=unique(.x[[diag]]), factor=nm, miRNA=m, test=NA, p=NA_real_)
          } else {
            res <- wilcox_or_ttest(.x[[m]], .x$F)
            data.frame(Diagnosis=unique(.x[[diag]]), factor=nm, miRNA=m, test=res$test, p=res$p)
          }
        }) %>% ungroup()

      tables[[paste0("Factor_", nm, "_", m)]] <- list(
        summary = summary_group(sub, m, c(diag, "F")),
        test = tab
      )
    }
  }

  # export all expression tables to one xlsx
  out_xlsx <- file.path(out_dir, "tables", "expression_suite_tables.xlsx")
  write_expression_workbook(tables, out_xlsx)

  list(tables = tables, plots = plots, out_xlsx = out_xlsx)
}

write_expression_workbook <- function(tables, path) {
  wb <- openxlsx::createWorkbook()
  for (nm in names(tables)) {
    sh <- str_sub(nm, 1, 31)
    openxlsx::addWorksheet(wb, sh)
    obj <- tables[[nm]]
    r <- 1
    for (piece in names(obj)) {
      openxlsx::writeData(wb, sh, x = data.frame(Section = piece), startRow = r, startCol = 1)
      r <- r + 1
      openxlsx::writeData(wb, sh, x = obj[[piece]], startRow = r, startCol = 1)
      r <- r + nrow(obj[[piece]]) + 3
    }
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
}