# prepare_10cancers_mean_deseq2_vst_inputs
.libPaths(c(.libPaths(), "D:/SIGLNCOS_revision_LiuMingyu/02_survival_reanalysis/envs/r_seurat44/lib/R/library"))
suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(hdf5r)
})

root <- "D:/SIGLNCOS_revision_LiuMingyu/02_survival_reanalysis"
gdc_dir <- file.path(root, "intermediate/gdc_release_45")
pair8_file <- file.path(root, "metadata/all_cancer_survival_inputs/legacy_survival_pairs_8cancers_with_compartment.csv")
pair2_file <- file.path(root, "results/tables/legacy_pair_audit.csv")
map_file <- file.path(root, "metadata/host_cell_to_music_compartment_mapping.csv")

out_matrix_dir <- file.path(root, "intermediate/deseq2_vst_10cancers_mean")
out_data_dir <- file.path(root, "intermediate/pair_score_lock_10cancers_mean_deseq2_vst")
out_meta_dir <- file.path(root, "metadata/pair_score_lock_10cancers_mean_deseq2_vst")
dir.create(out_matrix_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_meta_dir, recursive = TRUE, showWarnings = FALSE)

exploratory_cancers <- c("ICC", "UVM")
main_cancers <- c("BC", "BLCA", "CRC", "GBM", "HNSCC", "LUAD", "OV", "UCEC")

pairs8 <- fread(pair8_file)[
  ,
  .(
    cancer, project_id, legacy_host_cell_type, lncRNA_1, lncRNA_2,
    coindex_sign, music_compartment
  )
]

mapping <- fread(map_file)[, .(cancer, legacy_host_cell_type, music_compartment)]

pairs2 <- fread(pair2_file)[
  cancer %in% c("LUAD", "BLCA") & archived_survival_flag == 1,
  .(
    cancer,
    project_id = paste0("TCGA-", cancer),
    legacy_host_cell_type = host_cell_type,
    lncRNA_1,
    lncRNA_2,
    coindex_sign
  )
]
pairs2 <- merge(
  pairs2,
  mapping,
  by = c("cancer", "legacy_host_cell_type"),
  all.x = TRUE,
  all.y = FALSE
)
if (pairs2[is.na(music_compartment), .N] > 0) {
  stop("Missing LUAD/BLCA host-cell to MuSiC compartment mapping.")
}

pairs <- rbindlist(list(pairs8, pairs2), fill = TRUE)
setorder(pairs, cancer, legacy_host_cell_type, lncRNA_1, lncRNA_2)

fraction_dir_for <- function(cancer) {
  if (cancer %in% c("LUAD", "BLCA")) {
    file.path(root, "results/deconvolution/music_compartment_full")
  } else {
    file.path(root, "results/deconvolution/music_8cancers_compartment")
  }
}

pair_rows <- list()
pair_summary_rows <- list()
gene_rows <- list()
pair_idx <- 0L
sum_idx <- 0L
gene_idx <- 0L

for (ca in unique(pairs$cancer)) {
  cancer_pairs <- pairs[cancer == ca]
  project <- unique(cancer_pairs$project_id)
  if (length(project) != 1) stop("Expected one project per cancer for ", ca)
  project <- project[[1]]

  h5_path <- file.path(gdc_dir, project, "bulk_expression_primary_patient.h5")
  frac_path <- file.path(fraction_dir_for(ca), paste0(project, "_MuSiC_compartment_proportions.csv"))

  h5 <- H5File$new(h5_path, mode = "r")
  counts <- h5[["counts_unstranded"]][, ]
  genes <- h5[["gene_name"]][]
  patient_id <- h5[["patient_id"]][]
  sample_id <- h5[["sample_id"]][]
  h5$close_all()

  if (nrow(counts) == length(patient_id) && ncol(counts) == length(genes)) counts <- t(counts)
  stopifnot(nrow(counts) == length(genes), ncol(counts) == length(patient_id))
  rownames(counts) <- genes
  colnames(counts) <- patient_id

  targets <- sort(unique(c(cancer_pairs$lncRNA_1, cancer_pairs$lncRNA_2)))
  target_idx <- match(targets, genes)
  if (anyNA(target_idx)) {
    stop(sprintf(
      "Missing target genes in %s: %s",
      project,
      paste(targets[is.na(target_idx)], collapse = ", ")
    ))
  }

  keep <- rowSums(counts >= 10) >= 10
  dds <- DESeqDataSetFromMatrix(
    countData = round(counts[keep, , drop = FALSE]),
    colData = data.frame(one = rep(1, length(patient_id)), row.names = patient_id),
    design = ~1
  )
  dds <- estimateSizeFactors(dds)
  vsd <- vst(dds, blind = TRUE, nsub = min(1000, nrow(dds)))
  expr <- assay(vsd)[targets, , drop = FALSE]
  saveRDS(expr, file.path(out_matrix_dir, paste0(project, "_target_lncRNA_DESeq2_VST_matrix.rds")))

  frac <- fread(frac_path)
  cohort <- merge(
    data.table(patient_id = patient_id, sample_id = sample_id),
    frac,
    by = "patient_id",
    all.x = FALSE,
    all.y = FALSE
  )

  for (gene in targets) {
    x <- as.numeric(expr[gene, match(cohort$patient_id, colnames(expr))])
    sx <- sd(x)
    gene_idx <- gene_idx + 1L
    gene_rows[[gene_idx]] <- data.table(
      cancer = ca,
      project_id = project,
      gene = gene,
      patients = length(x),
      vst_mean = mean(x),
      vst_sd = sx,
      eligible_nonzero_variance = is.finite(sx) && sx > 0
    )
  }

  for (i in seq_len(nrow(cancer_pairs))) {
    rec <- cancer_pairs[i]
    g1 <- rec$lncRNA_1
    g2 <- rec$lncRNA_2
    x1 <- as.numeric(expr[g1, match(cohort$patient_id, colnames(expr))])
    x2 <- as.numeric(expr[g2, match(cohort$patient_id, colnames(expr))])
    sd1 <- sd(x1)
    sd2 <- sd(x2)

    if (!(is.finite(sd1) && sd1 > 0 && is.finite(sd2) && sd2 > 0)) {
      sum_idx <- sum_idx + 1L
      pair_summary_rows[[sum_idx]] <- data.table(
        cancer = ca,
        project_id = project,
        legacy_host_cell_type = rec$legacy_host_cell_type,
        music_compartment = rec$music_compartment,
        lncRNA_1 = g1,
        lncRNA_2 = g2,
        coindex_sign = rec$coindex_sign,
        pair_score_eligible = FALSE,
        reason = "single-gene zero variance"
      )
      next
    }

    if (!(rec$music_compartment %in% names(cohort))) {
      stop(sprintf("Compartment %s not found in fraction table for %s", rec$music_compartment, project))
    }

    z1 <- as.numeric(scale(x1))
    z2 <- as.numeric(scale(x2))
    raw <- (z1 + z2) / 2
    score_sd <- sd(raw)
    if (!(is.finite(score_sd) && score_sd > 0)) {
      sum_idx <- sum_idx + 1L
      pair_summary_rows[[sum_idx]] <- data.table(
        cancer = ca,
        project_id = project,
        legacy_host_cell_type = rec$legacy_host_cell_type,
        music_compartment = rec$music_compartment,
        lncRNA_1 = g1,
        lncRNA_2 = g2,
        coindex_sign = rec$coindex_sign,
        pair_score_eligible = FALSE,
        reason = "pair-score zero variance"
      )
      next
    }

    score <- as.numeric(scale(raw))
    fraction <- cohort[[rec$music_compartment]]
    pair_id <- paste(ca, rec$legacy_host_cell_type, g1, g2, sep = "|")

    pair_idx <- pair_idx + 1L
    pair_rows[[pair_idx]] <- data.table(
      pair_id = pair_id,
      cancer = ca,
      project_id = project,
      patient_id = cohort$patient_id,
      sample_id = cohort$sample_id,
      legacy_host_cell_type = rec$legacy_host_cell_type,
      music_compartment = rec$music_compartment,
      lncRNA_1 = g1,
      lncRNA_2 = g2,
      coindex_sign = rec$coindex_sign,
      normalization = "DESeq2_VST",
      pair_score_strategy = "unified_mean_all_pairs",
      pair_score_formula = "(z_gene1 + z_gene2) / 2",
      primary_model_type = ifelse(ca %in% exploratory_cancers, "pairscore_cellfraction", "multivariable"),
      gene1_vst = x1,
      gene2_vst = x2,
      gene1_z = z1,
      gene2_z = z2,
      pair_score_raw = raw,
      pair_score_z = score,
      pair_score_median_cutoff = median(score),
      host_cell_fraction = fraction,
      host_cell_fraction_per_0_10 = fraction / 0.10
    )

    sum_idx <- sum_idx + 1L
    pair_summary_rows[[sum_idx]] <- data.table(
      cancer = ca,
      project_id = project,
      legacy_host_cell_type = rec$legacy_host_cell_type,
      music_compartment = rec$music_compartment,
      lncRNA_1 = g1,
      lncRNA_2 = g2,
      coindex_sign = rec$coindex_sign,
      pair_score_eligible = TRUE,
      reason = "",
      primary_model_type = ifelse(ca %in% exploratory_cancers, "pairscore_cellfraction", "multivariable"),
      patients = nrow(cohort),
      pair_score_sd = score_sd,
      pair_score_median_cutoff = median(score),
      host_fraction_zero_fraction = mean(fraction == 0),
      host_fraction_median = median(fraction)
    )
  }

  message("Prepared ", ca, " (", project, "): ", uniqueN(cohort$patient_id), " patients, ", nrow(cancer_pairs), " pairs")
}

locked <- rbindlist(pair_rows, fill = TRUE)
gene_audit <- rbindlist(gene_rows, fill = TRUE)
pair_summary <- rbindlist(pair_summary_rows, fill = TRUE)

fwrite(locked, file.path(out_data_dir, "patient_pair_scores_and_host_fractions_10cancers_mean_deseq2_vst.csv.gz"))
fwrite(gene_audit, file.path(out_meta_dir, "pair_gene_expression_eligibility_10cancers_mean_deseq2_vst.csv"))
fwrite(pair_summary, file.path(out_meta_dir, "pair_score_lock_summary_10cancers_mean_deseq2_vst.csv"))
fwrite(
  data.table(
    cancer = c(main_cancers, exploratory_cancers),
    normalization = "DESeq2_VST",
    pair_score_strategy = "unified_mean_all_pairs",
    pair_score_formula = "(z_gene1 + z_gene2) / 2",
    primary_model_type = c(rep("multivariable", length(main_cancers)), rep("pairscore_cellfraction", length(exploratory_cancers)))
  ),
  file.path(out_meta_dir, "included_cancers_10cancers_mean_deseq2_vst.csv")
)

print(locked[, .(patients = uniqueN(patient_id), pairs = uniqueN(pair_id)), by = cancer][order(cancer)])


# run_survival_models_10cancers_mean_deseq2_vst
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

root <- "D:/SIGLNCOS_revision_LiuMingyu/02_survival_reanalysis"
score_file <- file.path(root, "intermediate/pair_score_lock_10cancers_mean_deseq2_vst/patient_pair_scores_and_host_fractions_10cancers_mean_deseq2_vst.csv.gz")
clinical_dir1 <- file.path(root, "intermediate/all_cancer_survival_inputs")
clinical_dir2 <- file.path(root, "metadata/gdc_release_45")
out_dir <- file.path(root, "results/survival_10cancers_mean_deseq2_vst")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scores <- fread(score_file)
project_map <- unique(scores[, .(cancer, project_id, primary_model_type)])

safe_cox <- function(formula, data) {
  tryCatch(
    coxph(formula, data = data, x = TRUE, y = TRUE, model = TRUE, ties = "efron"),
    error = function(e) e
  )
}

coefficient_table <- function(fit, pair_id, cancer, analysis, model_type) {
  if (inherits(fit, "error")) return(data.table())
  z <- summary(fit)$coefficients
  ci <- summary(fit)$conf.int
  data.table(
    pair_id = pair_id,
    cancer = cancer,
    analysis = analysis,
    model_type = model_type,
    term = rownames(z),
    beta = z[, "coef"],
    HR = ci[, "exp(coef)"],
    CI_low = ci[, "lower .95"],
    CI_high = ci[, "upper .95"],
    se = z[, "se(coef)"],
    z = z[, "z"],
    p = z[, "Pr(>|z|)"]
  )
}

model_audit <- function(fit, d, pair_id, cancer, model_type, rhs_base) {
  if (inherits(fit, "error")) {
    return(data.table(
      pair_id = pair_id,
      cancer = cancer,
      model_type = model_type,
      n = nrow(d),
      events = sum(d$OS_event),
      parameters = NA_integer_,
      events_per_parameter = NA_real_,
      AIC = NA_real_,
      concordance = NA_real_,
      model_covariates = rhs_base,
      error = conditionMessage(fit)
    ))
  }
  data.table(
    pair_id = pair_id,
    cancer = cancer,
    model_type = model_type,
    n = nrow(d),
    events = sum(d$OS_event),
    parameters = length(coef(fit)),
    events_per_parameter = sum(d$OS_event) / length(coef(fit)),
    AIC = AIC(fit),
    concordance = unname(summary(fit)$concordance[1]),
    model_covariates = rhs_base,
    error = NA_character_
  )
}

observed_levels_n <- function(x) length(unique(as.character(x[!is.na(x)])))

load_clinical <- function(project, cancer) {
  path1 <- file.path(clinical_dir1, paste0(project, "_clinical_expression_patient_master.csv"))
  path2 <- file.path(clinical_dir2, project, "clinical_expression_patient_master.csv")
  if (file.exists(path1)) {
    x <- fread(path1)
  } else if (file.exists(path2)) {
    x <- fread(path2)
  } else {
    stop("Missing clinical master for ", project)
  }
  if (!("cancer" %in% names(x))) x[, cancer := cancer]
  x[, age10 := age_years / 10]
  x[, sex_model := factor(tolower(sex), levels = c("female", "male"))]
  x[, stage_binary := fifelse(stage_group %in% c("I", "II"), "early",
                       fifelse(stage_group %in% c("III", "IV"), "advanced", NA_character_))]
  x[, stage_binary := factor(stage_binary, levels = c("early", "advanced"))]
  x
}

clinical_list <- list()
for (i in seq_len(nrow(project_map))) {
  cancer <- project_map$cancer[[i]]
  project <- project_map$project_id[[i]]
  clinical_list[[cancer]] <- load_clinical(project, cancer)
}

pair_ids <- unique(scores$pair_id)
coef_rows <- list()
audit_rows <- list()
cohort_rows <- list()
idx_out <- 0L

for (pid in pair_ids) {
  d <- copy(scores[pair_id == pid])
  cancer <- unique(d$cancer)
  model_type <- unique(d$primary_model_type)
  if (length(model_type) != 1) stop("Expected one primary model type per pair for ", pid)
  model_type <- model_type[[1]]

  d <- merge(
    d,
    clinical_list[[cancer]],
    by = c("patient_id", "sample_id", "cancer", "project_id"),
    all.x = TRUE
  )
  d <- d[
    is.finite(OS_days) & OS_days > 0 & OS_event %in% c(0, 1) &
      is.finite(pair_score_z) & is.finite(gene1_z) & is.finite(gene2_z) &
      is.finite(host_cell_fraction)
  ]
  setorder(d, patient_id)

  if (model_type == "multivariable") {
    covars <- c("age10")
    if (observed_levels_n(d$sex_model) >= 2) covars <- c(covars, "sex_model")
    if (observed_levels_n(d$stage_binary) >= 2) covars <- c(covars, "stage_binary")
    covars <- c(covars, "host_cell_fraction_per_0_10")
    needed <- c("OS_days", "OS_event", "pair_score_z", "gene1_z", "gene2_z", covars)
    d <- d[complete.cases(d[, ..needed])]
    analysis_pair <- "pair_full"
    analysis_g1 <- "gene1_full"
    analysis_g2 <- "gene2_full"
  } else {
    covars <- c("host_cell_fraction_per_0_10")
    needed <- c("OS_days", "OS_event", "pair_score_z", "gene1_z", "gene2_z", covars)
    d <- d[complete.cases(d[, ..needed])]
    analysis_pair <- "pair_reduced"
    analysis_g1 <- "gene1_reduced"
    analysis_g2 <- "gene2_reduced"
  }

  rhs_base <- paste(covars, collapse = " + ")

  if (nrow(d) == 0L) {
    audit_rows[[pid]] <- data.table(
      pair_id = pid,
      cancer = cancer,
      model_type = model_type,
      n = 0,
      events = 0,
      parameters = NA_integer_,
      events_per_parameter = NA_real_,
      AIC = NA_real_,
      concordance = NA_real_,
      model_covariates = rhs_base,
      error = "no complete cases"
    )
    next
  }

  fit_pair <- safe_cox(as.formula(paste("Surv(OS_days, OS_event) ~", rhs_base, "+ pair_score_z")), d)
  fit_g1 <- safe_cox(as.formula(paste("Surv(OS_days, OS_event) ~", rhs_base, "+ gene1_z")), d)
  fit_g2 <- safe_cox(as.formula(paste("Surv(OS_days, OS_event) ~", rhs_base, "+ gene2_z")), d)

  idx_out <- idx_out + 1L
  coef_rows[[idx_out]] <- coefficient_table(fit_pair, pid, cancer, analysis_pair, model_type)
  idx_out <- idx_out + 1L
  coef_rows[[idx_out]] <- coefficient_table(fit_g1, pid, cancer, analysis_g1, model_type)
  idx_out <- idx_out + 1L
  coef_rows[[idx_out]] <- coefficient_table(fit_g2, pid, cancer, analysis_g2, model_type)

  audit_rows[[pid]] <- model_audit(fit_pair, d, pid, cancer, model_type, rhs_base)
  cohort_rows[[pid]] <- unique(d[, .(
    pair_id = pid,
    cancer,
    project_id,
    normalization,
    pair_score_strategy,
    pair_score_formula,
    model_type,
    legacy_host_cell_type,
    music_compartment,
    lncRNA_1,
    lncRNA_2,
    coindex_sign,
    n = .N,
    events = sum(OS_event),
    model_covariates = rhs_base
  )])
}

coef_dt <- rbindlist(coef_rows, fill = TRUE)
audit_dt <- rbindlist(audit_rows, fill = TRUE)
cohort_dt <- rbindlist(cohort_rows, fill = TRUE)

pair_primary <- coef_dt[
  term == "pair_score_z" & analysis %in% c("pair_full", "pair_reduced")
]
pair_primary[, p_FDR_BH := p.adjust(p, method = "BH"), by = cancer]

cancer_summary <- merge(
  cohort_dt[, .(
    pairs = .N,
    median_n = as.numeric(median(n)),
    median_events = as.numeric(median(events)),
    model_type = unique(model_type)[1]
  ), by = cancer],
  audit_dt[, .(median_events_per_parameter = median(events_per_parameter, na.rm = TRUE)), by = cancer],
  by = "cancer",
  all = TRUE
)
cancer_summary[, pair_score_formula := "(z_gene1 + z_gene2) / 2"]
cancer_summary[, normalization := "DESeq2_VST"]

fwrite(coef_dt, file.path(out_dir, "all_coefficients_10cancers_mean_deseq2_vst.csv"))
fwrite(pair_primary, file.path(out_dir, "pair_primary_effects_with_FDR_10cancers_mean_deseq2_vst.csv"))
fwrite(audit_dt, file.path(out_dir, "pair_model_audit_10cancers_mean_deseq2_vst.csv"))
fwrite(cohort_dt, file.path(out_dir, "pair_complete_case_cohort_10cancers_mean_deseq2_vst.csv"))
fwrite(cancer_summary, file.path(out_dir, "cancer_model_scope_summary_10cancers_mean_deseq2_vst.csv"))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "R_sessionInfo.txt"))

print(pair_primary[, .(pairs = .N, nominal_p_lt_0_05 = sum(p < 0.05), FDR_lt_0_05 = sum(p_FDR_BH < 0.05)), by = .(cancer, model_type)][order(cancer)])


# build_single_column_forest_10cancers_mean_deseq2_vst
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

root <- "D:/SIGLNCOS_revision_LiuMingyu/02_survival_reanalysis"
surv_dir <- file.path(root, "results/survival_10cancers_mean_deseq2_vst")
out_dir <- file.path(root, "results/publication_package/figures/forest_10cancers_mean_deseq2_vst_single_column")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

coef <- fread(file.path(surv_dir, "all_coefficients_10cancers_mean_deseq2_vst.csv"))
eff <- fread(file.path(surv_dir, "pair_primary_effects_with_FDR_10cancers_mean_deseq2_vst.csv"))
cohort <- fread(file.path(surv_dir, "pair_complete_case_cohort_10cancers_mean_deseq2_vst.csv"))

compartment_order <- c("Immune", "Immune_hematopoietic", "Immune_myeloid", "Fibroblast", "Endothelial", "Epithelial_tumor", "Tumor", "Tumor_glial")
compartment_label <- c(
  Immune = "Immune",
  Immune_hematopoietic = "Immune / hematopoietic",
  Immune_myeloid = "Immune / myeloid",
  Fibroblast = "Fibroblast",
  Endothelial = "Endothelial",
  Epithelial_tumor = "Tumor / epithelial",
  Tumor = "Tumor",
  Tumor_glial = "Tumor / glial"
)

term_maps <- list(
  multivariable = c(
    pair_score_z = "co-lncRNA PairScore (mean z, per SD)",
    age10 = "Age (per 10 years)",
    sex_modelmale = "Sex (male vs female)",
    stage_binaryadvanced = "Stage (advanced vs early)",
    host_cell_fraction_per_0_10 = "Matched cell fraction (per 10%)"
  ),
  pairscore_cellfraction = c(
    pair_score_z = "co-lncRNA PairScore (mean z, per SD)",
    host_cell_fraction_per_0_10 = "Matched cell fraction (per 10%)"
  )
)

fmt_p <- function(x) ifelse(x < 1e-4, "<0.0001", formatC(x, format = "f", digits = 4))
fmt_n <- function(x) {
  y <- formatC(x, format = "f", digits = 2)
  sub("0+$", "", sub("\\.$", "", y))
}

sp <- tstrsplit(eff$pair_id, "|", fixed = TRUE)
eff[, `:=`(legacy_cell = sp[[2]], gene1 = sp[[3]], gene2 = sp[[4]])]
eff <- merge(
  eff,
  unique(cohort[, .(pair_id, model_type, music_compartment)]),
  by = c("pair_id", "model_type"),
  all.x = TRUE
)
selected <- eff[order(p, p_FDR_BH, gene1, gene2), head(.SD, 2), by = .(cancer, music_compartment)]

build_forest_source <- function() {
  rows <- list()
  idx <- 0L
  for (pid in selected$pair_id) {
    meta <- cohort[pair_id == pid][1]
    mt <- meta$model_type
    z <- coef[pair_id == pid & model_type == mt]
    tmap <- term_maps[[mt]]
    pair_analysis <- if (mt == "multivariable") "pair_full" else "pair_reduced"
    g1_analysis <- if (mt == "multivariable") "gene1_full" else "gene1_reduced"
    g2_analysis <- if (mt == "multivariable") "gene2_full" else "gene2_reduced"

    full <- z[
      analysis == pair_analysis & term %in% names(tmap),
      .(variable = tmap[term], HR, lo = CI_low, hi = CI_high, p, type = ifelse(term == "pair_score_z", "PairScore", "Covariate"))
    ]
    g1 <- z[analysis == g1_analysis & term == "gene1_z"][1]
    g2 <- z[analysis == g2_analysis & term == "gene2_z"][1]
    pair_name <- paste0("PairScore: ", meta$lncRNA_1, " | ", meta$lncRNA_2)
    full[type == "PairScore", variable := pair_name]

    rr <- rbind(
      full,
      data.table(variable = meta$lncRNA_1, HR = g1$HR, lo = g1$CI_low, hi = g1$CI_high, p = g1$p, type = "Single lncRNA"),
      data.table(variable = meta$lncRNA_2, HR = g2$HR, lo = g2$CI_low, hi = g2$CI_high, p = g2$p, type = "Single lncRNA"),
      fill = TRUE
    )
    ord <- c(pair_name, meta$lncRNA_1, meta$lncRNA_2, unname(tmap[names(tmap) != "pair_score_z"]))
    rr[, `:=`(
      pair_id = pid,
      cancer = meta$cancer,
      model_type = mt,
      compartment = meta$music_compartment,
      pair = paste(meta$lncRNA_1, meta$lncRNA_2, sep = " | "),
      n = meta$n,
      events = meta$events,
      order = match(variable, ord)
    )]
    idx <- idx + 1L
    rows[[idx]] <- rr
  }
  ans <- rbindlist(rows, fill = TRUE)
  ans[, compartment_display := fifelse(compartment %in% names(compartment_label), compartment_label[compartment], compartment)]
  ans
}

build_layout_table <- function(z) {
  comp_levels <- unique(z$compartment)
  comp_levels <- compartment_order[compartment_order %in% comp_levels]
  pair_order <- unique(z[, .(pair_id, compartment, pair)])
  pair_order[, compartment := factor(compartment, levels = comp_levels)]
  setorder(pair_order, compartment, pair)

  rows <- list()
  idx <- 0L
  pair_block <- 0L

  for (comp in comp_levels) {
    ids <- pair_order[as.character(compartment) == comp, pair_id]
    if (!length(ids)) next
    comp_label <- unique(z[compartment == comp, compartment_display])[[1]]

    idx <- idx + 1L
    rows[[idx]] <- data.table(
      row_type = "section",
      display_label = comp_label,
      pair_id = NA_character_,
      block_id = NA_integer_,
      text_color = "#222222",
      HR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_
    )

    for (pid in ids) {
      pair_block <- pair_block + 1L
      zz <- copy(z[pair_id == pid][order(order)])

      idx <- idx + 1L
      rows[[idx]] <- data.table(
        row_type = "pair_header",
        display_label = unique(zz$pair),
        pair_id = pid,
        block_id = pair_block,
        text_color = "#222222",
        HR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_
      )

      zz[, display_label := fifelse(type == "PairScore", sub("^PairScore: ", "", variable), variable)]
      zz[, block_id := pair_block]
      zz[, row_type := "data"]
      zz[, text_color := fifelse(type == "PairScore", "#B2182B", fifelse(type == "Single lncRNA", "#2166AC", "#222222"))]
      zz[, subgroup := fifelse(type == "PairScore", "co-lncRNA", fifelse(type == "Single lncRNA", "single lncRNA", "model factors"))]

      for (sg in c("co-lncRNA", "single lncRNA", "model factors")) {
        zzz <- zz[subgroup == sg]
        if (!nrow(zzz)) next
        idx <- idx + 1L
        rows[[idx]] <- data.table(
          row_type = "subsection",
          display_label = sg,
          pair_id = pid,
          block_id = pair_block,
          text_color = "#222222",
          HR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_
        )
        for (j in seq_len(nrow(zzz))) {
          idx <- idx + 1L
          rows[[idx]] <- zzz[j, .(row_type, display_label, pair_id, block_id, text_color, HR, lo, hi, p)]
        }
      }

      idx <- idx + 1L
      rows[[idx]] <- data.table(
        row_type = "spacer",
        display_label = "",
        pair_id = pid,
        block_id = pair_block,
        text_color = "#222222",
        HR = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_
      )
    }
  }

  d <- rbindlist(rows, fill = TRUE)
  step_map <- c(section = 1.55, pair_header = 1.10, subsection = 0.85, data = 0.86, spacer = 1.00)
  ypos <- numeric(nrow(d))
  cursor <- 0
  for (i in seq_len(nrow(d))) {
    cursor <- cursor + step_map[[d$row_type[i]]]
    ypos[i] <- cursor
  }
  d[, y := max(ypos) - ypos + 1]
  d[, p_text := ifelse(row_type == "data", fmt_p(p), "")]
  d[, hr_text := ifelse(row_type == "data", paste0(fmt_n(HR), " (", fmt_n(lo), "-", fmt_n(hi), ")"), "")]
  d
}

draw_single_column <- function(d, cancer, model_type) {
  x_min <- 0.5
  x_max <- 3.0
  x_breaks <- c(0.5, 1, 2, 3)
  shade_blocks <- unique(na.omit(d$block_id))
  shade_blocks <- shade_blocks[shade_blocks %% 2 == 0]
  shade_dt <- d[block_id %in% shade_blocks & row_type %in% c("pair_header", "subsection", "data")]
  shade_band <- shade_dt[, .(ymin = min(y) - 0.42, ymax = max(y) + 0.42), by = block_id]
  plot_dt <- copy(d[row_type == "data"])
  plot_dt[, `:=`(
    lo_plot = pmax(lo, x_min),
    hi_plot = pmin(hi, x_max),
    trunc_left = lo < x_min,
    trunc_right = hi > x_max,
    line_alpha = fifelse(text_color == "#222222", 0.60, 0.95),
    point_alpha = fifelse(text_color == "#222222", 0.75, 1.00)
  )]
  left_arrow_dt <- plot_dt[trunc_left == TRUE]
  right_arrow_dt <- plot_dt[trunc_right == TRUE]

  common <- theme_void() + theme(plot.margin = margin(4, 4, 4, 4))

  left <- ggplot() +
    geom_rect(data = shade_band, aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax), fill = "#EEF3F1", color = NA) +
    geom_text(data = d[row_type == "section"], aes(x = 0.02, y = y, label = display_label), hjust = 0, fontface = "bold", size = 4.1) +
    geom_text(data = d[row_type == "pair_header"], aes(x = 0.04, y = y, label = display_label), hjust = 0, fontface = "bold", size = 3.15) +
    geom_text(data = d[row_type == "subsection"], aes(x = 0.08, y = y, label = display_label), hjust = 0, fontface = "bold.italic", size = 2.95) +
    geom_text(data = d[row_type == "data"], aes(x = 0.12, y = y, label = display_label, colour = text_color), hjust = 0, size = 2.75) +
    scale_colour_identity() +
    coord_cartesian(xlim = c(0, 1), ylim = c(min(d$y) - 0.6, max(d$y) + 0.6), clip = "off") +
    common

  forest <- ggplot() +
    geom_rect(data = shade_band, aes(xmin = x_min, xmax = x_max, ymin = ymin, ymax = ymax), fill = "#EEF3F1", color = NA) +
    geom_vline(xintercept = 1, linetype = 2, colour = "#C85858", linewidth = 0.35) +
    geom_segment(data = plot_dt, aes(x = lo_plot, xend = hi_plot, y = y, yend = y, colour = text_color, alpha = line_alpha), linewidth = 0.42, lineend = "round") +
    geom_segment(data = left_arrow_dt, aes(x = x_min * 1.08, xend = x_min, y = y, yend = y, colour = text_color, alpha = line_alpha), linewidth = 0.42, lineend = "round", arrow = arrow(length = grid::unit(0.07, "in"), ends = "last", type = "closed")) +
    geom_segment(data = right_arrow_dt, aes(x = x_max / 1.08, xend = x_max, y = y, yend = y, colour = text_color, alpha = line_alpha), linewidth = 0.42, lineend = "round", arrow = arrow(length = grid::unit(0.07, "in"), ends = "last", type = "closed")) +
    geom_point(data = plot_dt, aes(x = pmin(pmax(HR, x_min), x_max), y = y, colour = text_color, alpha = point_alpha), shape = 15, size = 2.05) +
    scale_x_log10(breaks = x_breaks, labels = x_breaks) +
    scale_colour_identity() +
    scale_alpha_identity() +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(min(d$y) - 0.6, max(d$y) + 0.6), clip = "off") +
    labs(x = "HR (log scale)", y = NULL) +
    theme_classic(base_size = 8) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), plot.margin = margin(4, 4, 4, 4))

  right <- ggplot() +
    geom_rect(data = shade_band, aes(xmin = 0, xmax = 1, ymin = ymin, ymax = ymax), fill = "#EEF3F1", color = NA) +
    geom_text(data = d[row_type == "section"], aes(x = 0.05, y = y, label = "p.value"), hjust = 0, fontface = "bold", size = 3.0) +
    geom_text(data = d[row_type == "section"], aes(x = 0.42, y = y, label = "HR (95% CI)"), hjust = 0, fontface = "bold", size = 3.0) +
    geom_text(data = d[row_type == "data"], aes(x = 0.05, y = y, label = p_text, colour = text_color), hjust = 0, size = 2.6) +
    geom_text(data = d[row_type == "data"], aes(x = 0.42, y = y, label = hr_text, colour = text_color), hjust = 0, size = 2.6) +
    scale_colour_identity() +
    coord_cartesian(xlim = c(0, 1), ylim = c(min(d$y) - 0.6, max(d$y) + 0.6), clip = "off") +
    common

  subtitle_text <- if (model_type == "multivariable") {
    "Primary model: PairScore + age + sex + stage + matched cell fraction. PairScore = mean(z_gene1, z_gene2)."
  } else {
    "Exploratory primary model: PairScore + matched cell fraction only. PairScore = mean(z_gene1, z_gene2)."
  }

  left + forest + right + plot_layout(widths = c(1.55, 1.15, 1.0)) +
    plot_annotation(title = paste0(cancer, " survival forest"), subtitle = subtitle_text)
}

src <- build_forest_source()
fwrite(src, file.path(out_dir, "integrated_forest_source_10cancers_mean_deseq2_vst.csv"))

summary_rows <- list()
idx <- 0L
for (ca in unique(src$cancer)) {
  z <- copy(src[cancer == ca])
  mt <- unique(z$model_type)[[1]]
  d <- build_layout_table(z)
  fig <- draw_single_column(d, ca, mt)
  stem <- paste0("Result5_", ca, "_DESeq2_VST_mean_single_column_forest")
  h <- max(8.5, nrow(d) * 0.23)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), fig, width = 10.8, height = h, units = "in")
  ggsave(file.path(out_dir, paste0(stem, ".png")), fig, width = 10.8, height = h, units = "in", dpi = 300)
  ggsave(file.path(out_dir, paste0(stem, ".tiff")), fig, width = 10.8, height = h, units = "in", dpi = 450, compression = "lzw")
  fwrite(d, file.path(out_dir, paste0(stem, "_displayed_rows.csv")))

  idx <- idx + 1L
  summary_rows[[idx]] <- unique(z[, .(
    cancer,
    model_type,
    displayed_pairs = uniqueN(pair_id),
    displayed_compartments = uniqueN(compartment),
    pair_score_formula = "(z_gene1 + z_gene2) / 2"
  )])
}

fwrite(rbindlist(summary_rows, fill = TRUE), file.path(out_dir, "forest_display_summary_10cancers_mean_deseq2_vst.csv"))



# build_10cancers_mean_km

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(patchwork)
})

root <- "D:/SIGLNCOS_revision_LiuMingyu/02_survival_reanalysis"
score_file <- file.path(root, "intermediate/pair_score_lock_10cancers_mean_deseq2_vst/patient_pair_scores_and_host_fractions_10cancers_mean_deseq2_vst.csv.gz")
forest_dir <- file.path(root, "results/publication_package/figures/forest_10cancers_mean_deseq2_vst_single_column")
out_dir <- file.path(root, "results/publication_package/figures/km_10cancers_mean_deseq2_vst")
clinical_dir1 <- file.path(root, "intermediate/all_cancer_survival_inputs")
clinical_dir2 <- file.path(root, "metadata/gdc_release_45")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

scores <- fread(score_file)
cohort_meta <- fread(file.path(root, "results/survival_10cancers_mean_deseq2_vst/pair_complete_case_cohort_10cancers_mean_deseq2_vst.csv"))

load_clinical <- function(project, cancer) {
  path1 <- file.path(clinical_dir1, paste0(project, "_clinical_expression_patient_master.csv"))
  path2 <- file.path(clinical_dir2, project, "clinical_expression_patient_master.csv")
  if (file.exists(path1)) {
    x <- fread(path1)
  } else if (file.exists(path2)) {
    x <- fread(path2)
  } else {
    stop("Missing clinical master for ", project)
  }
  if (!("cancer" %in% names(x))) x[, cancer := cancer]
  x[, age10 := age_years / 10]
  x[, sex_model := factor(tolower(sex), levels = c("female", "male"))]
  x[, stage_binary := fifelse(stage_group %in% c("I", "II"), "early",
                       fifelse(stage_group %in% c("III", "IV"), "advanced", NA_character_))]
  x[, stage_binary := factor(stage_binary, levels = c("early", "advanced"))]
  x
}

find_cutpoint <- function(d, minprop = 0.20) {
  x <- sort(unique(d$pair_score_z))
  lo <- as.numeric(quantile(d$pair_score_z, minprop, na.rm = TRUE))
  hi <- as.numeric(quantile(d$pair_score_z, 1 - minprop, na.rm = TRUE))
  cand <- x[x >= lo & x <= hi]
  if (!length(cand)) return(median(d$pair_score_z, na.rm = TRUE))
  chisq <- sapply(cand, function(cc) {
    grp <- d$pair_score_z > cc
    if (mean(grp) < minprop || mean(!grp) < minprop) return(NA_real_)
    fit <- tryCatch(survdiff(Surv(OS_days, OS_event) ~ grp, d), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    unname(fit$chisq)
  })
  if (all(!is.finite(chisq))) return(median(d$pair_score_z, na.rm = TRUE))
  cand[[which.max(chisq)]]
}

get_displayed_pairs <- function() {
  files <- list.files(forest_dir, pattern = "_displayed_rows\\.csv$", full.names = TRUE)
  rows <- rbindlist(lapply(files, fread), fill = TRUE)
  unique(rows[!is.na(pair_id) & pair_id != "", .(pair_id)])
}

km_plot <- function(d, pair_label, pval, x_breaks) {
  fit <- survfit(Surv(OS_days, OS_event) ~ group, data = d)
  sm <- summary(fit)
  curve_dt <- data.table(
    time = sm$time,
    surv = sm$surv,
    strata = as.character(sm$strata)
  )
  curve_dt[, group := sub("^group=", "", strata)]
  curve_dt <- rbindlist(list(
    data.table(time = 0, surv = 1, group = "Low"),
    data.table(time = 0, surv = 1, group = "High"),
    curve_dt[, .(time, surv, group)]
  ), fill = TRUE)

  risk_sm <- summary(fit, times = x_breaks, extend = TRUE)
  risk_dt <- data.table(
    time = risk_sm$time,
    n_risk = risk_sm$n.risk,
    strata = as.character(risk_sm$strata)
  )
  risk_dt[, group := sub("^group=", "", strata)]
  risk_dt[, y := fifelse(group == "Low", 2, 1)]
  risk_header <- data.table(time = x_breaks, n_risk = as.character(x_breaks), group = "Header", y = 3)
  risk_lbl <- data.table(time = min(x_breaks), n_risk = c("No. at risk", "Low", "High"), group = c("HeaderLbl", "LowLbl", "HighLbl"), y = c(3, 2, 1))

  pal <- c(Low = "#2E86DE", High = "#E74C3C")
  p1 <- ggplot(curve_dt, aes(x = time, y = surv, colour = group)) +
    geom_step(linewidth = 0.9) +
    scale_colour_manual(values = pal) +
    scale_y_continuous(limits = c(0, 1), expand = c(0.02, 0.02)) +
    scale_x_continuous(breaks = x_breaks, expand = c(0.01, 0.01)) +
    labs(
      title = pair_label,
      subtitle = paste0("Log-rank p = ", ifelse(pval < 1e-4, "<0.0001", formatC(pval, format = "f", digits = 4))),
      x = "Time (days)",
      y = "Overall survival probability",
      colour = ""
    ) +
    theme_classic(base_size = 10) +
    theme(legend.position = "top")

  p2 <- ggplot() +
    geom_text(data = risk_header, aes(x = time, y = y, label = n_risk), size = 3.2, fontface = "bold") +
    geom_text(data = risk_dt, aes(x = time, y = y, label = n_risk, colour = group), size = 3.2) +
    geom_text(data = risk_lbl, aes(x = time, y = y, label = n_risk), hjust = 0, size = 3.2, fontface = ifelse(risk_lbl$group == "HeaderLbl", "bold", "plain")) +
    scale_colour_manual(values = pal, guide = "none") +
    scale_y_continuous(limits = c(0.5, 3.5), breaks = NULL) +
    scale_x_continuous(breaks = x_breaks, labels = rep("", length(x_breaks)), expand = c(0.05, 0.05)) +
    theme_void() +
    theme(plot.margin = margin(0, 10, 5, 45))

  p1 / p2 + plot_layout(heights = c(3.3, 1.1))
}

displayed_pairs <- get_displayed_pairs()$pair_id
summary_rows <- list(); group_rows <- list(); idx <- 0L; gidx <- 0L
plot_list_by_cancer <- list()

for (pid in displayed_pairs) {
  meta <- cohort_meta[pair_id == pid][1]
  cancer <- meta$cancer
  project <- meta$project_id
  clinical <- load_clinical(project, cancer)
  d <- copy(scores[pair_id == pid])
  d <- merge(d, clinical, by = c("patient_id", "sample_id", "cancer", "project_id"), all.x = TRUE)

  covars <- trimws(unlist(strsplit(meta$model_covariates, "\\+")))
  covars <- covars[nzchar(covars)]
  needed <- unique(c("OS_days", "OS_event", "pair_score_z", "host_cell_fraction", covars))
  d <- d[is.finite(OS_days) & OS_days > 0 & OS_event %in% c(0, 1) & complete.cases(d[, ..needed])]
  if (nrow(d) < 10) next

  cut <- find_cutpoint(d, minprop = 0.20)
  d[, group := factor(ifelse(pair_score_z > cut, "High", "Low"), levels = c("Low", "High"))]
  lr <- survdiff(Surv(OS_days, OS_event) ~ group, data = d)
  pval <- 1 - pchisq(unname(lr$chisq), df = 1)
  xmax <- max(d$OS_days, na.rm = TRUE)
  x_breaks <- pretty(c(0, xmax), n = 5)
  x_breaks <- x_breaks[x_breaks >= 0 & x_breaks <= xmax]
  if (length(x_breaks) < 3) x_breaks <- unique(round(seq(0, xmax, length.out = 5)))
  if (min(x_breaks) > 0) x_breaks <- c(0, x_breaks)

  pair_label <- paste0(meta$lncRNA_1, " | ", meta$lncRNA_2)
  fig <- km_plot(d, pair_label, pval, x_breaks)
  stem <- paste0("KM_", cancer, "_", gsub("[|]", "_", pid))
  ggsave(file.path(out_dir, paste0(stem, ".png")), fig, width = 6.6, height = 6.3, dpi = 300)
  ggsave(file.path(out_dir, paste0(stem, ".pdf")), fig, width = 6.6, height = 6.3)

  plot_list_by_cancer[[cancer]] <- c(plot_list_by_cancer[[cancer]], list(fig))

  idx <- idx + 1L
  summary_rows[[idx]] <- d[, .(
    pair_id = pid, cancer = cancer, project_id = project, model_type = meta$model_type,
    lncRNA_1 = meta$lncRNA_1, lncRNA_2 = meta$lncRNA_2,
    cutoff = cut,
    n_total = .N,
    events_total = sum(OS_event),
    n_low = sum(group == "Low"), n_high = sum(group == "High"),
    events_low = sum(OS_event[group == "Low"]), events_high = sum(OS_event[group == "High"]),
    logrank_p = pval
  )][1]

  gidx <- gidx + 1L
  group_rows[[gidx]] <- d[, .(pair_id = pid, patient_id, pair_score_z, cutoff = cut, group, OS_days, OS_event)]
}

for (ca in names(plot_list_by_cancer)) {
  combo <- wrap_plots(plot_list_by_cancer[[ca]], ncol = 2)
  ggsave(file.path(out_dir, paste0("KM_", ca, "_combined_mean_deseq2_vst.pdf")), combo, width = 13.2, height = max(6.6, 6.6 * ceiling(length(plot_list_by_cancer[[ca]]) / 2)))
}

summary_dt <- rbindlist(summary_rows, fill = TRUE)
groups_dt <- rbindlist(group_rows, fill = TRUE)

fwrite(summary_dt, file.path(out_dir, "KM_cutpoint_summary_10cancers_mean_deseq2_vst.csv"))
fwrite(groups_dt, file.path(out_dir, "KM_group_assignments_10cancers_mean_deseq2_vst.csv"))

print(summary_dt[, .(pairs = .N), by = cancer][order(cancer)])
