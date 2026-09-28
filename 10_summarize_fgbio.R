# Credentials--------------------------------------------------------------

#
# Automated analysis of liquid biopsy data from fgbio workflow (HiLiB)
# Author: Axel Künstner
# CHIP-code based on implementation by Leopold Schawe
#

# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Axel Künstner
# Part of the LION panel analysis code: https://github.com/kunstner/UKSH_liquidbiopsy
# For research use only; not validated as a medical device.

# Description:
#   Per-sample summarization of cfDNA sequencing data processed through the
#   fgbio UMI-based error correction pipeline. For each sample the script
#   produces an annotated variant table (VarDict and Mutect2), applies
#   Panel-of-Normals and noise filtering, classifies CHIP variants (M-CHIP,
#   L-CHIP pathogenic, L-CHIP putative), estimates surrogate tumor fraction
#   (TF) following Husain et al. (2022, JCO Precis Oncol), and exports
#   results as Excel workbooks alongside QC metrics and visualizations.
#
# Run:
#   Interactive : source('10_summarize_fgbio.R') or run in RStudio; uses BATCHID defined below
#   Command line: Rscript 10_summarize_fgbio.R <batch_id>  (overrides BATCHID)
#
# Dependencies:
#   00_helper.R, 01_chip_classifier.R
#
# Output:
#   <BATCH>/fgbio_tables/Results_<sample>_Table.xlsx   per-sample variant tables
#   <BATCH>/fgbio_tables/QC_Summary_<BATCH>.xlsx       batch-level QC summary
#                                                       (sheet TF_Summary added)
#   <BATCH>/fgbio_tables/QC_Summary_<BATCH>_plots.pdf
#   <BATCH>/fgbio_tables/QC_Advanced_<BATCH>_plots.pdf
#   <BATCH>/fgbio_tables/Results_<sample>_heatmap.pdf
#   <BATCH>/fgbio_tables/Results_<sample>_vaf.pdf
#
# batch execution (batches 1 - NUMBER_BATCHES)

# for i in $(seq -w 1 NUMBER_BATCHES); do
# Rscript 10_summarize_fgbio.R "batch_a0${i}"
# done

# Libraries ---------------------------------------------------------------

pacman::p_load(tidyverse, openxlsx2, patchwork)

# Load helper functions
source('00_helper.R')

# functions for CHiP classification
source('01_chip_classifier.R')

# Setup -------------------------------------------------------------------

BATCHID <- 'batch_a001'
PATHARCHIVE <- '../'

# BATCH can be passed as a command-line argument:
#   Rscript 10_summarize_fgbio.R batch_a011
# If no argument is given, the value defined below is used as fallback.

args <- commandArgs(trailingOnly = TRUE)

BATCH  <- if (length(args) >= 1) args[1] else BATCHID

WRKFLW <- 'fgbio'
dir    <- paste0('../', BATCH, '/', WRKFLW, '_results')

cat("Processing batch:", BATCH, "\n")

DBPATH_MCHIP      <- 'DB/M-CHIP_Genes_modified_v1.2.tsv'
DBPATH_LCHIP_PATH <- 'DB/L-CHIP_Genes_pathogenic_modified_v1.2.tsv'
DBPATH_LCHIP_PUT  <- 'DB/L-CHIP_Genes_putative_modified_v1.2.tsv'

options("openxlsx2.na" = "")

# Gene Panel, Databases ---------------------------------------------------

# CHiP Classification
chip_db_mchip  <- load_chip_db(DBPATH_MCHIP,      "M-CHIP")
chip_db_lpath  <- load_chip_db(DBPATH_LCHIP_PATH, "L-CHIP_pathogenic")
chip_db_lput   <- load_chip_db(DBPATH_LCHIP_PUT,  "L-CHIP_putative")

# Panel genes
gene_list <- read_delim(file = '../data/Gene_list.txt', 
                        show_col_types = FALSE) |> 
    dplyr::mutate(Gene = case_when(
        Gene == 'FAM46C (TENT5C)' ~ 'TENT5C',
        .default = Gene
    )) |> 
    arrange(Gene)

# Data DrugBank
approved_drugs <- read_delim(file = 'EXTERNAL_DB/DrugBank_5.1.12/all drug links.csv', 
                             show_col_types = FALSE, name_repair = make.names)
protein_drug <- read_delim(file = 'EXTERNAL_DB/DrugBank_5.1.12/drugbank_all_target_polypeptide_ids/pharmacologically_active.csv', 
                           show_col_types = FALSE, name_repair = make.names)

# Panel of Normals data
pon_background <- openxlsx2::read_xlsx(file = 'Panel_of_normal_fgbio.xlsx', 
                                       sheet = 'Background') |> 
    dplyr::mutate(IDX = paste(Chromosome, Start_Position, Ref, Alt, sep = '_'))
pon_noise  <- openxlsx2::read_xlsx(file = 'Panel_of_normal_fgbio.xlsx', 
                                   sheet = 'Noise') |> 
    dplyr::mutate(IDX = paste(Chromosome, Start_Position, Ref, Alt, sep = '_'))

# Data File Lists ---------------------------------------------------------

# Input files
ff_mutect2 <- list.files(path = dir, pattern = "*mutect2.MAF.gz$", 
                         full.names = TRUE, recursive = TRUE)
ff_vardict <- list.files(path = dir, pattern = "*VarDict.MAF.gz$", 
                         full.names = TRUE, recursive = TRUE)

# QC files
ff_qc_01 <- list.files(path = dir, pattern = "^01_pcr_metrics.txt.gz$", 
                       full.names = TRUE, recursive = TRUE)
ff_qc_02 <- list.files(path = dir, pattern = "^03_consensus_pcr_metrics.txt.gz$", 
                       full.names = TRUE, recursive = TRUE)
ff_qc_03 <- list.files(path = dir, pattern = "^04_pcr_metrics.txt.gz$", 
                       full.names = TRUE, recursive = TRUE)
ff_qc_04 <- list.files(path = dir, pattern = "^04_stats.insert.txt.gz$", 
                       full.names = TRUE, recursive = TRUE)
ff_fusi <- list.files(path = dir, pattern = "^fusions_factera.txt.gz$", 
                      full.names = TRUE, recursive = TRUE)

# Create outdir
outdir <- gsub('results', 'tables', dir)
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Configuration -----------------------------------------------------------

# Variant classes
v_classes <- c('Frame_Shift_Del', 'Frame_Shift_Ins',
               'In_Frame_Del', 'In_Frame_Ins',
               'Missense_Mutation', 'Nonsense_Mutation', 'Nonstop_Mutation',
               'Splice_Region', 'Splice_Site')

# VAF thresholds for filtering
vaf_thresholds <- list(
    "0.0050" = 0.005,
    "0.0025" = 0.0025,
    "0.0005" = 0.0005
)

# QC processing configuration
qc_configs <- list(
    list(files = ff_qc_01, name = "Mapping", insert_stats = FALSE),
    list(files = ff_qc_02, name = "Error correction", insert_stats = FALSE),
    list(files = ff_qc_03, name = "Remapping clipped reads", insert_stats = FALSE),
    list(files = ff_qc_04, name = "Remapping coverage", insert_stats = TRUE)
)

# Data Completeness Validation -------------------------------------------

# Check that all file lists have the same length
file_lists <- list(
    "VarDict" = ff_vardict,
    "Mutect2" = ff_mutect2,
    "QC_01" = ff_qc_01,
    "QC_02" = ff_qc_02,
    "QC_03" = ff_qc_03,
    "QC_04" = ff_qc_04,
    "Fusions" = ff_fusi
)

n_samples <- length(ff_qc_01)
cat("Expected number of samples:", n_samples, "\n")

# Check file count consistency
for(file_type in names(file_lists)) {
    if(length(file_lists[[file_type]]) != n_samples) {
        stop(paste0("ERROR: Inconsistent file count for ", file_type, 
                    ". Expected: ", n_samples, 
                    ", Found: ", length(file_lists[[file_type]])))
    }
}

# Check that all required files exist
missing_files_summary <- data.frame(
    Sample_Index = integer(),
    Sample_Name = character(),
    Missing_Files = character(),
    stringsAsFactors = FALSE
)

for(i in 1:n_samples) {
    sample <- gsub('fullrun_|/vep.VarDict.MAF.gz', '', ff_vardict[i])
    sample <- tail(str_split(sample, pattern = '/')[[1]], 1)
    
    missing_files <- c()
    for(file_type in names(file_lists)) {
        if(!file.exists(file_lists[[file_type]][i])) {
            missing_files <- c(missing_files, file_type)
        }
    }
    
    if(length(missing_files) > 0) {
        missing_files_summary <- rbind(missing_files_summary, 
                                       data.frame(Sample_Index = i,
                                                  Sample_Name = sample,
                                                  Missing_Files = paste(missing_files, collapse = ", ")))
    }
}

# Report and stop if critical files missing
if(nrow(missing_files_summary) > 0) {
    cat("ERROR: Missing required files detected:\n")
    print(missing_files_summary)
    stop("Pipeline cannot proceed with missing files. Please check your data directory.")
} else {
    cat(" All required files present for", n_samples, "samples\n")
}

# Collector for batch-level TF summary ------------------------------------

# Initialise an empty list; one tibble row per sample, bound after the loop.
all_tf_results <- vector("list", length(ff_mutect2))

# Main Processing Loop ----------------------------------------------------

for(i in 1:length(ff_mutect2)) {
    
    mylist <- list()
    
    sample <- gsub('fullrun_|/vep.VarDict.MAF.gz', '', ff_vardict[i])
    sample <- tail(str_split(sample, pattern = '/')[[1]], 1)
    
    print(paste0('Screening sample ', sample, " (", i, ")"))
    
    # Process QC Results --------------------------------------------------
    
    qc_results <- map2_dfr(qc_configs, 
                           map(qc_configs, ~.x$files[i]),
                           ~read_qc_file(.y, .x$name, .x$insert_stats))
    
    mylist[['QC']] <- qc_results
    
    # Process Variant Files -----------------------------------------------
    
    # VarDict variants
    vardict_variants <- process_variants(
        file_path = ff_vardict[i],
        sample_name = sample,
        gene_list = gene_list,
        pon_background = pon_background,
        pon_noise = pon_noise,
        v_classes = v_classes,
        protein_drug = protein_drug,
        approved_drugs = approved_drugs,
        is_mutect = FALSE
    )
    
    mylist[['VarDict Variants']] <- vardict_variants
    
    # VarDict filtered variants at different thresholds
    for(threshold_name in names(vaf_thresholds)) {
        sheet_name <- paste('VarDict Variants Filt.', threshold_name)
        mylist[[sheet_name]] <- apply_variant_filters(
            variants = vardict_variants,
            vaf_threshold = vaf_thresholds[[threshold_name]],
            alt_count_threshold = 2,
            max_af_threshold = 0.001,
            filter_pon = TRUE,
            filter_variants = TRUE,
            v_classes = v_classes
        )
    }
    
    # VarDict CHIP
    mylist[['VarDict CHIP M-CHIP']] <- 
        classify_chip(variants = vardict_variants, chip_db = chip_db_mchip,
                      vaf_threshold = 0.01) |> 
        dplyr::select(-CHIP_criteria, -Drug_ID, -Drug)
    mylist[['VarDict CHIP L-CHIP pathogenic']] <- 
        classify_chip(variants = vardict_variants, chip_db = chip_db_lpath,
                      vaf_threshold = 0.01) |> 
        dplyr::select(-CHIP_criteria, -Drug_ID, -Drug)
    mylist[['VarDict CHIP L-CHIP putative']] <- 
        classify_chip(variants = vardict_variants, chip_db = chip_db_lput,
                      vaf_threshold = 0.01) |> 
        dplyr::select(-CHIP_criteria, -Drug_ID, -Drug)
    
    # Mutect2 variants
    mutect2_variants <- process_variants(
        file_path = ff_mutect2[i],
        sample_name = sample,
        gene_list = gene_list,
        pon_background = pon_background,
        pon_noise = pon_noise,
        v_classes = v_classes,
        protein_drug = protein_drug,
        approved_drugs = approved_drugs,
        is_mutect = TRUE
    )
    
    mylist[['Mutect2 Variants']] <- mutect2_variants
    
    # Mutect2 PASS variants
    mylist[['Mutect2 Variants PASS']] <- apply_variant_filters(
        variants = mutect2_variants,
        max_af_threshold = 0.001,
        filter_pon = TRUE,
        filter_variants = TRUE,
        filter_pass = TRUE,
        v_classes = v_classes
    )
    
    # Mutect2 CHIP
    mylist[['Mutect2 CHIP M-CHIP']] <- 
        classify_chip(variants = mutect2_variants, chip_db = chip_db_mchip,
                      vaf_threshold = 0.01) |> 
        dplyr::select(-CHIP_criteria, -Drug_ID, -Drug)
    mylist[['Mutect2 CHIP L-CHIP pathogenic']] <- 
        classify_chip(variants = mutect2_variants, chip_db = chip_db_lpath,
                      vaf_threshold = 0.01) |> 
        dplyr::select(-CHIP_criteria, -Drug_ID, -Drug)
    mylist[['Mutect2 CHIP L-CHIP putative']] <- 
        classify_chip(variants = mutect2_variants, chip_db = chip_db_lput,
                      vaf_threshold = 0.01) |> 
        dplyr::select(-CHIP_criteria, -Drug_ID, -Drug)
    
    # Fusions
    mylist[['Fusions']] <- read_delim(file = ff_fusi[i], show_col_types = FALSE)
    
    # Tumor Fraction Estimation -------------------------------------------
    #
    # Surrogate TF following Husain et al. (2022): max somatic VAF after
    # excluding CHIP variants. Uses the 0.0050 VarDict filtered table as
    # input (PoN-filtered, noise-filtered, v_classes restricted) to keep
    # the somatic pool clean. Mutect2 PASS concordance is checked for the
    # max-VAF variant as an additional confidence flag.
    #
    # A KDE-based TF estimate (experimental) is appended when enough
    # variants are available (see estimate_tf_kde() in 00_helper.R).
    
    tf_input <- mylist[['VarDict Variants Filt. 0.0050']] |>
        dplyr::filter(is.na(t_vaf_above_noise) | t_vaf_above_noise == TRUE)
    
    tf_row <- estimate_tumor_fraction(
        vardict_filtered = tf_input,
        # vardict_filtered = mylist[['VarDict Variants Filt. 0.0050']],
        mutect2_pass     = mylist[['Mutect2 Variants PASS']],
        chip_mchip       = mylist[['VarDict CHIP M-CHIP']],
        chip_lpath       = mylist[['VarDict CHIP L-CHIP pathogenic']],
        chip_lput        = mylist[['VarDict CHIP L-CHIP putative']],
        sample_name      = sample
    )
    
    # Reconstruct somatic pool for KDE ----------------------------------------
    # The somatic pool is not exposed by estimate_tumor_fraction() to keep its
    # return value a clean one-row tibble; it is reconstructed here using the
    # same CHIP exclusion and germline guard logic so both functions operate on
    # identical input.
    chip_idx_loop <- dplyr::bind_rows(
        mylist[['VarDict CHIP M-CHIP']],
        mylist[['VarDict CHIP L-CHIP pathogenic']],
        mylist[['VarDict CHIP L-CHIP putative']]
    ) |>
        dplyr::mutate(idx = paste(Chromosome, Start_Position,
                                  Tumor_Seq_Allele1, Tumor_Seq_Allele2, sep = "_")) |>
        dplyr::pull(idx) |>
        unique()
    
    somatic_loop <- tf_input |>
        dplyr::mutate(
            .idx = paste(Chromosome, Start_Position,
                         Tumor_Seq_Allele1, Tumor_Seq_Allele2, sep = "_"),
            t_AF = suppressWarnings(as.numeric(t_AF))
        ) |>
        dplyr::filter(!.idx %in% chip_idx_loop,
                      !is.na(t_AF),
                      t_AF < 0.35)
    
    tf_row <- dplyr::bind_cols(tf_row, estimate_tf_kde(somatic_loop))
    
    mylist[['Tumor Fraction']] <- tf_row
    
    # Store for batch-level summary
    all_tf_results[[i]] <- tf_row
    
    cat(sprintf(
        "  TF max: %s  |  p80: %s  |  mean: %s (sd: %s)  |  median: %s  |  kde: %s  |  tier: %s  |  M2: %s\n",
        ifelse(is.na(tf_row$TF_max),    "NA", sprintf("%.2f%%", tf_row$TF_max    * 100)),
        ifelse(is.na(tf_row$TF_p80),    "NA", sprintf("%.2f%%", tf_row$TF_p80    * 100)),
        ifelse(is.na(tf_row$TF_mean),   "NA", sprintf("%.2f%%", tf_row$TF_mean   * 100)),
        ifelse(is.na(tf_row$TF_sd),     "NA", sprintf("%.2f%%", tf_row$TF_sd     * 100)),
        ifelse(is.na(tf_row$TF_median), "NA", sprintf("%.2f%%", tf_row$TF_median * 100)),
        ifelse(is.na(tf_row$TF_kde),    "NA", sprintf("%.2f%%", tf_row$TF_kde    * 100)),
        tf_row$TF_tier,
        ifelse(is.na(tf_row$TF_mutect2_concordant), "NA",
               as.character(tf_row$TF_mutect2_concordant))
    ))
    
    # Export Results ------------------------------------------------------
    
    mylist <- mylist[c(
        "QC",
        "Tumor Fraction",
        "VarDict Variants",
        "VarDict Variants Filt. 0.0050",
        "VarDict Variants Filt. 0.0025",
        "VarDict Variants Filt. 0.0005",
        "VarDict CHIP M-CHIP",
        "VarDict CHIP L-CHIP pathogenic",
        "VarDict CHIP L-CHIP putative",
        "Mutect2 Variants",
        "Mutect2 Variants PASS",
        "Mutect2 CHIP M-CHIP",
        "Mutect2 CHIP L-CHIP pathogenic",
        "Mutect2 CHIP L-CHIP putative",
        "Fusions"
    )]
    
    openxlsx2::write_xlsx(x = mylist, 
                          file = paste0(outdir, '/Results_', sample, '_Table.xlsx'), 
                          first_row = TRUE)
    
    # Create Visualizations -----------------------------------------------
    
    # Heatmap
    create_gene_heatmap(
        vardict_variants = mylist[['VarDict Variants Filt. 0.0050']],
        mutect2_variants = mylist[['Mutect2 Variants PASS']],
        gene_list = gene_list,
        sample_name = sample,
        outdir = outdir
    )
    
    # VAF plots
    tryCatch(
        expr = {
            create_vaf_plots(
                variants   = vardict_variants,
                sample_name = sample,
                outdir     = outdir,
                binwidth   = 0.005
            )
        },
        error = function(e) {
            warning("create_vaf_plots() failed for sample ", sample, ": ", conditionMessage(e))
        }
    )
}

# QC Summary Across All Samples ------------------------------------------

cat("\nGenerating QC summary across all samples...\n")

# Collect all QC results
all_qc_results <- data.frame()

for(i in 1:length(ff_qc_01)) {
    sample <- gsub('fullrun_|/vep.VarDict.MAF.gz', '', ff_vardict[i])
    sample <- tail(str_split(sample, pattern = '/')[[1]], 1)
    
    # Read QC data for this sample
    sample_qc <- map2_dfr(qc_configs, 
                          map(qc_configs, ~.x$files[i]),
                          ~read_qc_file(.y, .x$name, .x$insert_stats)) |> 
        dplyr::mutate(Sample = sample)
    
    all_qc_results <- dplyr::bind_rows(all_qc_results, sample_qc)
}

# Create summary statistics
qc_summary <- all_qc_results |> 
    dplyr::group_by(Processing) |> 
    dplyr::summarise(
        n_samples = dplyr::n(),
        # Total reads
        mean_total_reads = round(mean(TOTAL_READS, na.rm = TRUE), 0),
        median_total_reads = round(median(TOTAL_READS, na.rm = TRUE), 0),
        min_total_reads = min(TOTAL_READS, na.rm = TRUE),
        max_total_reads = max(TOTAL_READS, na.rm = TRUE),
        # Coverage metrics
        mean_coverage = round(mean(MEAN_TARGET_COVERAGE, na.rm = TRUE), 0),
        median_coverage = round(median(MEAN_TARGET_COVERAGE, na.rm = TRUE), 0),
        min_coverage = min(MEAN_TARGET_COVERAGE, na.rm = TRUE),
        max_coverage = max(MEAN_TARGET_COVERAGE, na.rm = TRUE),
        # Insert size metrics  
        mean_insert_size = round(mean(MEAN_INSERT_SIZE, na.rm = TRUE), 0),
        median_insert_size = round(median(MEAN_INSERT_SIZE, na.rm = TRUE), 0),
        .groups = 'drop'
    )

# Create detailed sample-wise QC table
qc_detailed <- all_qc_results |> 
    dplyr::select(Sample, Processing, TOTAL_READS, MEAN_TARGET_COVERAGE, 
                  MEDIAN_TARGET_COVERAGE, MEAN_INSERT_SIZE, STANDARD_DEVIATION) |> 
    tidyr::pivot_wider(names_from = Processing, 
                       values_from = c(TOTAL_READS, MEAN_TARGET_COVERAGE, 
                                       MEDIAN_TARGET_COVERAGE, MEAN_INSERT_SIZE, 
                                       STANDARD_DEVIATION),
                       names_sep = "_")

# Bind per-sample TF rows into a batch-level summary table ----------------
#
# TF_tier distribution is printed to console and saved as a dedicated sheet
# in the QC workbook. Samples with TF >= 10% and an empty somatic variant
# table are flagged as "confident negative" (Husain et al. 2022).

tf_summary <- dplyr::bind_rows(all_tf_results) |>
    dplyr::mutate(
        TF_max_pct    = dplyr::if_else(!is.na(TF_max),    round(TF_max    * 100, 2), NA_real_),
        TF_p80_pct    = dplyr::if_else(!is.na(TF_p80),    round(TF_p80    * 100, 2), NA_real_),
        TF_mean_pct   = dplyr::if_else(!is.na(TF_mean),   round(TF_mean   * 100, 2), NA_real_),
        TF_sd_pct     = dplyr::if_else(!is.na(TF_sd),     round(TF_sd     * 100, 2), NA_real_),
        TF_median_pct = dplyr::if_else(!is.na(TF_median), round(TF_median * 100, 2), NA_real_),
        # Flag samples where the max-VAF variant is a likely outlier relative
        # to the bulk VAF distribution. A ratio > 3 indicates that TF_max is
        # driven by a single high-VAF variant (e.g. a dominant CHIP clone or
        # a single clonal driver) rather than reflecting the overall somatic
        # VAF level. In these cases the TF tier should be interpreted with
        # caution and the TF_max_gene inspected manually.
        TF_max_outlier = dplyr::case_when(
            is.na(TF_max) | is.na(TF_p80) ~ NA,
            TF_p80 == 0                   ~ NA,
            TF_max / TF_p80 > 3           ~ TRUE,
            .default                      = FALSE
        )
    ) |>
    dplyr::select(
        Sample,
        TF_max, TF_max_pct,
        TF_p80, TF_p80_pct,
        TF_mean, TF_mean_pct,
        TF_sd, TF_sd_pct,
        TF_median, TF_median_pct,
        TF_kde, TF_kde_floor, TF_kde_n_above,
        TF_tier, TF_max_outlier,
        TF_n_somatic, TF_max_gene, TF_max_hgvsp,
        TF_mutect2_concordant
    )

# Print tier distribution to console
cat("\n=== TUMOR FRACTION TIER DISTRIBUTION ===\n")
print(dplyr::count(tf_summary, TF_tier, name = "n_samples"))

# Export QC summary
qc_export_list <- list(
    "QC_Summary"  = qc_summary,
    "QC_Detailed" = qc_detailed,
    "QC_All_Raw"  = all_qc_results,
    "TF_Summary"  = tf_summary
)

openxlsx2::write_xlsx(x = qc_export_list, 
                      file = paste0(outdir, '/QC_Summary_', BATCH, '.xlsx'), 
                      first_row = TRUE)

# Print summary to console
cat("\n=== QC SUMMARY ACROSS", length(ff_mutect2), "SAMPLES ===\n")
print(qc_summary)

# Check for potential issues
cat("\n=== QUALITY CHECKS ===\n")

# Low coverage samples
low_coverage_samples <- all_qc_results |> 
    dplyr::filter(Processing == "Mapping", MEAN_TARGET_COVERAGE < 100) |> 
    dplyr::pull(Sample)

if(length(low_coverage_samples) > 0) {
    cat("  Samples with low coverage (<100x):", paste(low_coverage_samples, collapse = ", "), "\n")
} else {
    cat(" All samples have adequate coverage (>=100x)\n")
}

# Low read count samples  
low_read_samples <- all_qc_results |> 
    dplyr::filter(Processing == "Mapping", TOTAL_READS < 1000000) |> 
    dplyr::pull(Sample)

if(length(low_read_samples) > 0) {
    cat(" Samples with low read count (<1M):", paste(low_read_samples, collapse = ", "), "\n")
} else {
    cat(" All samples have adequate read count (>=1M)\n")
}

# TF tier summary
cat("\n=== TF TIER SUMMARY ===\n")
undetectable_n <- sum(tf_summary$TF_tier == "undetectable")
low_tf_n       <- sum(tf_summary$TF_tier == "<1%")
mid_tf_n       <- sum(tf_summary$TF_tier == "1%-10%")
high_tf_n      <- sum(tf_summary$TF_tier == ">=10%")
cat(sprintf(
    " Undetectable: %d  |  <1%%: %d  |  1%%-10%%: %d  |  >=10%%: %d\n",
    undetectable_n, low_tf_n, mid_tf_n, high_tf_n
))
if(high_tf_n > 0) {
    high_tf_samples <- tf_summary |>
        dplyr::filter(TF_tier == ">=10%") |>
        dplyr::pull(Sample)
    cat(" High-TF samples (>=10%, near 100% sensitivity):",
        paste(high_tf_samples, collapse = ", "), "\n")
}
if(undetectable_n > 0) {
    cat(" NOTE:", undetectable_n,
        "sample(s) with undetectable TF — negative results are uninformative",
        "without tissue confirmation.\n")
}

cat("\n QC Summary saved to:", paste0(outdir, '/QC_Summary_', BATCH, '.xlsx\n'))

# Generate QC Plots -------------------------------------------------------

cat("Generating QC plots...\n")

qc_plots <- create_qc_plots(all_qc_results, BATCH, outdir)

cat(" QC plots saved to:", paste0(outdir, '/QC_Summary_', BATCH, '_plots.pdf\n'))

cat("Generating advanced QC plots...\n")

# Panel 1: On-Target Efficiency
# 
# (ON_TARGET_BASES / PF_BASES) x 100
# Shows percentage of sequencing effort hitting targets
# Higher is better for panel efficiency
# 
# Panel 2: Coverage Uniformity
# 
# MEDIAN_TARGET_COVERAGE / MEAN_TARGET_COVERAGE
# Values closer to 1.0 = more uniform coverage
# Red dashed line at 1.0 for reference
# 
# Panel 3: Read Retention
# 
# (current_reads / initial_reads) x 100
# Shows data loss through processing steps
# Tracks efficiency of error correction/remapping
# 
# Panel 4: Coverage Variability
# 
# log10(MAX_TARGET_COVERAGE - MIN_TARGET_COVERAGE + 1)
# Lower values = more consistent coverage across targets
# Log scale to handle wide ranges

advanced_qc_plots <- create_advanced_qc_plots(all_qc_results, BATCH, outdir)

cat(" Advanced QC plots saved to:", paste0(outdir, '/QC_Advanced_', BATCH, '_plots.pdf\n'))

# Copy results to archive folder ------------------------------------------

archive_dir <- file.path(PATHARCHIVE, 'liq_bio_runs', BATCH)

if (!dir.exists(archive_dir)) dir.create(archive_dir, recursive = TRUE)

file.copy(
    from      = list.files(outdir, full.names = TRUE),
    to        = archive_dir,
    overwrite = TRUE
)
