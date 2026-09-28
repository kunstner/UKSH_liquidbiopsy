# Credentials--------------------------------------------------------------

#
# Helper Functions for cfDNA Analysis Pipeline
# Author: Axel Künstner

# Function ACMG Classification --------------------------------------------

# Classification function taken from MIRACUM pipe
acmg <- function(df) {
    
    df <- data.frame(CLIN_SIG = df)
    
    # Split by ,
    df$CLIN_SIG <- str_split(df$CLIN_SIG, ",")
    
    # Remove multi annotations
    df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "[a-z_]+/[a-z_]+", ""))
    df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "^pathogenic$", "5"))
    df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "^likely_pathogenic$", "4"))
    df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "^conflicting_interpretations_of_pathogenicity$", "4*"))
    df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "^uncertain_significance$", "3"))
    df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "^likely_benign$", "2"))
    df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace_all(x, "^benign$", "1"))

    # Remove annotations without acmg relevance
    df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_remove(x, "^other$"))
    df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_remove(x, "^drug_response$"))
    df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_remove(x, "^risk_factor$"))
    df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_remove(x, "^not_provided$"))
    
    df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace(x, "^$", "0"))
    prio <- lapply(df$CLIN_SIG, function(x) str_replace(x, "4\\*", "4"))
    argmax <- unlist(lapply(prio, function(x) which.max(as.numeric(x))))
    df$CLIN_SIG <- lapply(df$CLIN_SIG, function(x) str_replace(x, "0", "."))
    
    return(mapply(function(x,y) x[y], df$CLIN_SIG, argmax))
}

# Drug Interactions Function ---------------------------------------------

drug_interactions <- function(vars, v_classes, protein_drug, approved_drugs) {
    vars <- vars |> 
        dplyr::mutate(Drug_ID = '', Drug = '')
    
    gens <- vars |> 
        dplyr::filter(Variant_Classification %in% v_classes) |> 
        dplyr::pull(Hugo_Symbol) |> 
        unique()
    
    for(ith in gens) {
        drugIDs <- protein_drug |> 
            dplyr::filter(Gene.Name %in% ith) |> 
            dplyr::pull(Drug.IDs) |> 
            str_split(., pattern = '; ') |> 
            unlist()
        
        drugNames <- approved_drugs |> 
            dplyr::filter(DrugBank.ID %in% drugIDs) |> 
            dplyr::pull(Name)
        
        vars$Drug_ID[vars$Hugo_Symbol == ith] <- str_flatten_comma(string = drugIDs)
        vars$Drug[vars$Hugo_Symbol == ith] <- str_flatten_comma(string = drugNames)
    }
    return(vars)
}


# Read QC Files Function -------------------------------------------------

read_qc_file <- function(file_path, processing_name, include_insert_stats = FALSE) {
    qc_data <- read_delim(file = file_path, 
                          skip = 6, n_max = 1, show_col_types = FALSE) |> 
        dplyr::mutate(Processing = processing_name)
    
    if (!include_insert_stats) {
        # Standard QC files with coverage metrics
        qc_data <- qc_data |> 
            dplyr::mutate(MEAN_TARGET_COVERAGE = round(MEAN_TARGET_COVERAGE, 0),
                          MEAN_INSERT_SIZE = NA, 
                          STANDARD_DEVIATION = NA)
    } else {
        # Insert stats file - different structure
        qc_data <- qc_data |> 
            dplyr::mutate(MEAN_INSERT_SIZE = round(MEAN_INSERT_SIZE, 0),
                          STANDARD_DEVIATION = round(STANDARD_DEVIATION, 0),
                          PF_BASES = NA, 
                          ON_TARGET_BASES = NA, 
                          TOTAL_READS = NA,
                          MEAN_TARGET_COVERAGE = NA, 
                          MEDIAN_TARGET_COVERAGE = NA, 
                          MIN_TARGET_COVERAGE = NA, 
                          MAX_TARGET_COVERAGE = NA)
    }
    
    return(qc_data |> 
               dplyr::select(Processing, TOTAL_READS, PF_BASES, ON_TARGET_BASES, 
                             MEAN_TARGET_COVERAGE, MEDIAN_TARGET_COVERAGE, 
                             MIN_TARGET_COVERAGE, MAX_TARGET_COVERAGE,
                             MEAN_INSERT_SIZE, STANDARD_DEVIATION))
}

# HGVSp Normalization ------------------------------------------------------

#' Normalize HGVSp-style protein notation before any downstream matching
#'
#' Two independent, idempotent fixes:
#'   1. Decode the URL-encoded equals sign VEP sometimes emits ("%3D" -> "=")
#'   2. Collapse an explicit same-residue substitution ("p.T125T" or
#'      "p.Thr125Thr") into correct HGVS synonymous notation ("p.T125=")
#'      -- covers annotations that wrote ref==alt letters instead of "=".
#'
#' @param x Character vector (HGVSc, HGVSp, or HGVSp_Short values)
#' @return Character vector, same length, normalized
.normalize_hgvsp <- function(x) {
    x <- stringr::str_replace_all(x, "%3D", "=")
    # 1-letter code: p.T125T -> p.T125=
    x <- stringr::str_replace(x, "^(p\\.)([A-Z])([0-9]+)\\2$", "\\1\\2\\3=")
    # 3-letter code: p.Thr125Thr -> p.Thr125=
    x <- stringr::str_replace(x, "^(p\\.)([A-Z][a-z]{2})([0-9]+)\\2$", "\\1\\2\\3=")
    x
}

# Process Variant File Function ------------------------------------------

process_variants <- function(file_path, sample_name, gene_list, pon_background, pon_noise,
                             v_classes, protein_drug, approved_drugs, is_mutect = FALSE) {
    
    variants <- read_delim(file = file_path, skip = 1, show_col_types = FALSE, na = '#N/A',
                           name_repair = if(!is_mutect) 'unique_quiet' else 'minimal') |> 
        dplyr::filter(Hugo_Symbol %in% gene_list$Gene) |>
        dplyr::mutate(Tumor_Sample_Barcode = sample_name)
    
    # Handle VarDict specific columns
    if (!is_mutect) {
        variants <- variants |> 
            dplyr::mutate(t_AF = t_alt_count/t_depth) |> 
            dplyr::select(-ends_with('77'), -ends_with('136'))
    }
    
    # Normalize HGVS notation: decode "%3D" and collapse same-AA
    # substitutions into "p.###=" notation; strip stray whitespace
    variants <- variants |>
        dplyr::mutate(dplyr::across(
            dplyr::any_of(c("HGVSc", "HGVSp", "HGVSp_Short")),
            \(x) x |> .normalize_hgvsp() |> stringr::str_squish()
        ))
    
    
    # Common processing
    variants <- variants |> 
        dplyr::mutate(ACMG_Class = acmg(CLIN_SIG)) |>
        dplyr::mutate(IDX = paste(Chromosome, Start_Position, 
                                  Tumor_Seq_Allele1, 
                                  Tumor_Seq_Allele2, sep = '_')) |> 
        dplyr::mutate(PoN_back = IDX %in% pon_background$IDX) |> 
        dplyr::left_join(pon_noise |> 
                             dplyr::select(IDX, vaf_noise, n_samples), 
                         by = 'IDX') |> 
        dplyr::mutate(t_vaf_above_noise = case_when(
            as.numeric(t_AF) > vaf_noise ~ TRUE,
            is.na(vaf_noise) ~ NA,
            .default = FALSE)) |> 
        dplyr::select(-IDX)
    
    # Add drug interactions
    variants <- drug_interactions(variants, v_classes, protein_drug, approved_drugs)
    
    return(variants)
}

# Apply Variant Filters Function ----------------------------------------

apply_variant_filters <- function(variants, vaf_threshold = NULL, alt_count_threshold = 2,
                                  max_af_threshold = 0.001, filter_pon = TRUE, 
                                  filter_variants = TRUE, v_classes = NULL,
                                  filter_pass = FALSE, ch_genes = NULL) {
    
    filtered <- variants
    
    # VAF filter
    if (!is.null(vaf_threshold)) {
        filtered <- filtered |> dplyr::filter(t_AF > vaf_threshold)
    }
    
    # Alt count filter
    if (!is.null(alt_count_threshold)) {
        filtered <- filtered |> dplyr::filter(t_alt_count > alt_count_threshold)
    }
    
    # Population frequency filter
    if (!is.null(max_af_threshold)) {
        filtered <- filtered |> dplyr::filter(is.na(MAX_AF) | MAX_AF <= max_af_threshold)
    }
    
    # PoN filter
    if (filter_pon) {
        filtered <- filtered |> dplyr::filter(PoN_back == FALSE)
    }
    
    # PASS filter (for Mutect2)
    if (filter_pass) {
        filtered <- filtered |> dplyr::filter(FILTER == 'PASS')
    }
    
    # Variant classification filter
    if (filter_variants && !is.null(v_classes)) {
        filtered <- filtered |> dplyr::filter(Variant_Classification %in% v_classes)
    }
    
    # CH genes filter
    if (!is.null(ch_genes)) {
        filtered <- filtered |> dplyr::filter(Hugo_Symbol %in% ch_genes$Gene)
    }
    
    return(filtered)
}

# Create Gene Heatmap Function -------------------------------------------

create_gene_heatmap <- function(vardict_variants, mutect2_variants, gene_list, sample_name, outdir) {
    
    tmp1 <- vardict_variants |>
        dplyr::filter(PoN_back == FALSE) |>
        dplyr::filter(is.na(t_vaf_above_noise) | t_vaf_above_noise == TRUE) |>
        dplyr::pull(Hugo_Symbol)
    
    tmp2 <- mutect2_variants |>
        dplyr::filter(PoN_back == FALSE) |>
        dplyr::filter(is.na(t_vaf_above_noise) | t_vaf_above_noise == TRUE) |>
        dplyr::pull(Hugo_Symbol)
    
    genes_long <- gene_list |>
        dplyr::select(-Inclusion) |>
        dplyr::mutate(
            Mutated = dplyr::if_else(Gene %in% c(tmp1, tmp2), "Mutated", "WT"),
            VarDict = dplyr::if_else(Gene %in% tmp1,          "Mutated", "WT"),
            Mutect2 = dplyr::if_else(Gene %in% tmp2,          "Mutated", "WT")) |>
        tidyr::pivot_longer(
            cols      = c(Mutated, VarDict, Mutect2),
            names_to  = "Caller",
            values_to = "Status") |>
        dplyr::mutate(
            Gene = factor(Gene, levels = rev(sort(unique(Gene)))),
            # Numeric x position with a manual gap of 0.6 between Mutated and callers
            x_pos = dplyr::case_when(
                Caller == "Mutated" ~ 1,
                Caller == "VarDict" ~ 2.6,
                Caller == "Mutect2" ~ 3.6))
    
    # Fixed cell size in inches
    cell_size   <- 0.15
    n_genes     <- dplyr::n_distinct(genes_long$Gene)
    n_cols      <- 3L
    plot_width  <- n_cols * cell_size + 1.15
    plot_height <- n_genes * cell_size + 0.7
    
    hmap_p <- ggplot(genes_long,
                     aes(x = x_pos, y = Gene, fill = interaction(Status, Caller))) +
        geom_tile(color = "white", linewidth = 0.3, width = 0.9, height = 0.9) +
        scale_fill_manual(
            values = c(
                "WT.Mutated"      = "grey90",
                "Mutated.Mutated" = "#B22222",
                "WT.VarDict"      = "grey90",
                "Mutated.VarDict" = "#4DBBD5",
                "WT.Mutect2"      = "grey90",
                "Mutated.Mutect2" = "#3C5488")) +
        scale_x_continuous(
            breaks = c(1, 2.6, 3.6),
            labels = c("Mutated", "VarDict", "Mutect2"),
            expand = expansion(add = 0.6)) +
        coord_fixed(ratio = 1) +
        labs(title = sample_name, x = NULL, y = NULL) +
        theme_classic() +
        theme(
            legend.position = "none",
            axis.text.x     = element_text(size = 8, angle = 45, hjust = 1),
            axis.text.y     = element_text(size = 8, face = "italic"),
            axis.line       = element_blank(),
            axis.ticks      = element_blank(),
            plot.title      = element_text(size = 9, hjust = 0.5))
    
    ggsave(
        filename = paste0(outdir, '/Results_', sample_name, '_heatmap.pdf'),
        height   = plot_height,
        width    = plot_width + 1,
        plot     = hmap_p)
    
    #    return(hmap_p)
}

# Create VAF Plots Function ----------------------------------------------

create_vaf_plots <- function(variants, sample_name, outdir, binwidth = 0.005) {
    
    hist_base <- variants |> 
        ggplot(aes(x = t_AF)) + 
        geom_histogram(binwidth = binwidth, color = "black", fill = "grey75") +
        geom_density(alpha = 0.2, fill = "#FF6666") +
        geom_vline(xintercept = c(0.01, 0.05), 
                   color = 'cornflowerblue', lty = 'dashed') +
        ggtitle(sample_name) + 
        xlab('Variant Allele Frequency') +
        theme_classic()
    
    hist_combined <- hist_base + xlim(0, 0.05) + 
        hist_base + xlim(0, 0.25) + 
        hist_base + 
        hist_base + xlim(0.9, 1) +
        plot_layout(ncol = 4)
    
    ggsave(filename = paste0(outdir, '/Results_', sample_name, '_vaf.pdf'), 
           height = 3, width = 12, plot = hist_combined)
    
    return(hist_combined)
}

# Create QC Summary Plots Function ---------------------------------------

create_qc_plots <- function(all_qc_results, batch_name, outdir) {
    
    # Filter for processing stages with coverage data
    coverage_stages <- c("Mapping", "Error correction", "Remapping clipped reads")
    
    coverage_data <- all_qc_results |> 
        dplyr::filter(Processing %in% coverage_stages) |> 
        dplyr::filter(!is.na(MEAN_TARGET_COVERAGE)) |> 
        dplyr::mutate(Processing = factor(Processing, levels = coverage_stages))
    
    # Create connected violin/boxplot
    p1 <- coverage_data |> 
        ggplot(aes(x = Processing, y = MEAN_TARGET_COVERAGE)) +
        geom_violin(fill = "lightblue", alpha = 0.6, trim = TRUE) +
        geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.alpha = 0) +
        geom_line(aes(group = Sample), alpha = 0.3, color = "darkgray") +
        geom_point(alpha = 0.6, size = 1.5, color = "darkblue") +
        scale_y_log10(labels = scales::label_number(accuracy = 1)) +
        labs(title = paste("Mean Target Coverage"),
             y = "Mean Target Coverage (x) [log scale]",
             x = "") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 12, face = "bold"))
    
    # Create total reads plot
    reads_data <- all_qc_results |> 
        dplyr::filter(Processing %in% coverage_stages) |> 
        dplyr::filter(!is.na(TOTAL_READS)) |> 
        dplyr::mutate(Processing = factor(Processing, levels = coverage_stages),
                      TOTAL_READS_M = TOTAL_READS / 1e6)  # Convert to millions
    
    p2 <- reads_data |> 
        ggplot(aes(x = Processing, y = TOTAL_READS_M)) +
        geom_violin(fill = "lightgreen", alpha = 0.6, trim = TRUE) +
        geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.alpha = 0) +
        geom_line(aes(group = Sample), alpha = 0.3, color = "darkgray") +
        geom_point(alpha = 0.6, size = 1.5, color = "darkgreen") +
        scale_y_log10() +
        labs(title = paste("Total Reads"),
             y = "Total Reads (Millions) [log scale]",
             x = "") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 12, face = "bold"))
    
    # Create insert size plot (only for remapping coverage stage)
    insert_data <- all_qc_results |> 
        dplyr::filter(Processing == "Remapping coverage") |> 
        dplyr::filter(!is.na(MEAN_INSERT_SIZE))
    
    if(nrow(insert_data) > 0) {
        p3 <- insert_data |> 
            ggplot(aes(x = "", y = MEAN_INSERT_SIZE)) +
            geom_violin(fill = "salmon", alpha = 0.6, trim = TRUE) +
            geom_boxplot(width = 0.3, fill = "white", alpha = 0.8) +
            geom_jitter(alpha = 0.6, size = 2, color = "darkred", width = 0.1) +
            labs(title = paste("Insert Sizes"),
                 y = "Mean Insert Size (bp)",
                 x = "") +
            theme_classic() +
            theme(plot.title = element_text(size = 12, face = "bold"))
        
        # Combine plots
        layout <- '
        AABBC
        '
        combined_plot <- p2 + p1 + p3 + 
            plot_layout(design = layout) + 
            plot_annotation( paste("Processed data - ", BATCH),
                             theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
    } else {
        # Combine plots without insert size
        combined_plot <- (p1 / p2)
    }
    
    # Save combined plot
    ggsave(filename = paste0(outdir, '/QC_Summary_', batch_name, '_plots.pdf'), 
           height = 6, width = 9, plot = combined_plot)
    
    return(list(summary = combined_plot))
}

# Create Advanced QC Plots Function --------------------------------------
create_advanced_qc_plots <- function(all_qc_results, batch_name, outdir) {
    
    # Filter for processing stages with coverage data
    coverage_stages <- c("Mapping", "Error correction", "Remapping clipped reads")
    
    # Calculate derived metrics
    qc_enhanced <- all_qc_results |> 
        dplyr::filter(Processing %in% coverage_stages) |> 
        dplyr::filter(!is.na(MEAN_TARGET_COVERAGE), !is.na(ON_TARGET_BASES), !is.na(PF_BASES)) |> 
        dplyr::mutate(Processing = factor(Processing, levels = coverage_stages)) |> 
        dplyr::mutate(
            # On-target efficiency
            On_Target_Efficiency = (ON_TARGET_BASES / PF_BASES) * 100,
            # Coverage uniformity
            Coverage_Uniformity = MEDIAN_TARGET_COVERAGE / MEAN_TARGET_COVERAGE,
            # Coverage range (log scale to handle large values)
            Coverage_Range_Log = log10(MAX_TARGET_COVERAGE - MIN_TARGET_COVERAGE + 1)
        ) |> 
        # Calculate read retention by sample
        dplyr::group_by(Sample) |> 
        dplyr::mutate(
            Initial_Reads = dplyr::first(TOTAL_READS),
            Read_Retention = (TOTAL_READS / Initial_Reads) * 100
        ) |> 
        dplyr::ungroup()
    
    # Plot 1: On-target efficiency
    p1 <- qc_enhanced |> 
        ggplot(aes(x = Processing, y = On_Target_Efficiency)) +
        geom_violin(fill = "lightcoral", alpha = 0.6, trim = TRUE) +
        geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.alpha = 0) +
        geom_line(aes(group = Sample), alpha = 0.3, color = "darkgray") +
        geom_point(alpha = 0.6, size = 1.5, color = "darkred") +
        labs(title = "On-Target Efficiency",
             y = "On-Target Efficiency (%)",
             x = "") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 10, face = "bold"))
    
    # Plot 2: Coverage uniformity
    p2 <- qc_enhanced |> 
        ggplot(aes(x = Processing, y = Coverage_Uniformity)) +
        geom_violin(fill = "lightgoldenrod", alpha = 0.6, trim = TRUE) +
        geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.alpha = 0) +
        geom_line(aes(group = Sample), alpha = 0.3, color = "darkgray") +
        geom_point(alpha = 0.6, size = 1.5, color = "darkgoldenrod") +
        geom_hline(yintercept = 1.0, linetype = "dashed", color = "red", alpha = 0.7) +
        labs(title = "Coverage Uniformity",
             y = "Median/Mean Coverage Ratio",
             x = "") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 10, face = "bold"))
    
    # Plot 3: Read retention
    p3 <- qc_enhanced |> 
        ggplot(aes(x = Processing, y = Read_Retention)) +
        geom_violin(fill = "lightseagreen", alpha = 0.6, trim = TRUE) +
        geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.alpha = 0) +
        geom_line(aes(group = Sample), alpha = 0.3, color = "darkgray") +
        geom_point(alpha = 0.6, size = 1.5, color = "darkseagreen") +
        labs(title = "Read Retention",
             y = "Read Retention (%)",
             x = "") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 10, face = "bold"))
    
    # Plot 4: Coverage range (variability)
    p4 <- qc_enhanced |> 
        ggplot(aes(x = Processing, y = Coverage_Range_Log)) +
        geom_violin(fill = "mediumpurple", alpha = 0.6, trim = TRUE) +
        geom_boxplot(width = 0.2, fill = "white", alpha = 0.8, outlier.alpha = 0) +
        geom_line(aes(group = Sample), alpha = 0.3, color = "darkgray") +
        geom_point(alpha = 0.6, size = 1.5, color = "darkslateblue") +
        labs(title = "Coverage Variability",
             y = "Coverage Range [log10(max-min+1)]",
             x = "") +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 10, face = "bold"))
    
    # Combine all advanced plots
    layout <- '
    ABCD
    '
    advanced_combined <- p1 + p2 + p3 + p4 +
        plot_layout(design = layout) + 
        plot_annotation(title = paste("Advanced QC Metrics -", batch_name),
                        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
    
    # Add overall title
    advanced_final <- advanced_combined + 
        plot_annotation(title = paste("Advanced QC Metrics -", batch_name),
                        theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5)))
    
    # Save advanced plot
    ggsave(filename = paste0(outdir, '/QC_Advanced_', batch_name, '_plots.pdf'), 
           height = 6, width = 12, plot = advanced_final)
    
    return(advanced_final)
}

# Estimate Tumor Fraction Function ---------------------------------------

#' Estimate surrogate tumor fraction (TF) from somatic variant VAFs
#'
#' Implements the variant-based TF estimation strategy of Husain et al. (2022,
#' JCO Precis Oncol): the highest allele fraction among non-germline,
#' non-CHIP somatic variants is used as a surrogate for ctDNA tumor fraction.
#' CHIP variants (from all three chip_db tables) are excluded before the
#' maximum is taken, mirroring the exclusion logic described in the paper.
#'
#' Additional descriptive statistics of the somatic VAF distribution are
#' reported alongside the max to characterise the spread of variant allele
#' frequencies and support interpretation of the TF estimate:
#'   - TF_p80    : 80th percentile VAF — a robust upper bound less sensitive
#'                 to single outlier variants than the max
#'   - TF_mean   : arithmetic mean VAF
#'   - TF_sd     : standard deviation of VAFs; a large SD relative to the mean
#'                 indicates a heterogeneous VAF distribution where the max
#'                 may not represent the dominant clone
#'   - TF_median : median VAF; consistently near the noise floor in low-TF
#'                 samples, useful as a noise floor indicator
#'
#' A large gap between TF_max and TF_p80 (or TF_mean/TF_median) suggests the
#' max-VAF variant is an outlier and should be interpreted cautiously.
#'
#' The function also checks whether the max-VAF variant is independently
#' supported by Mutect2 PASS calls (concordance flag), which increases
#' confidence in the estimate.
#'
#' TF tiers follow the paper's analytical thresholds:
#'   "undetectable"  : no somatic variants after CHIP and germline exclusion
#'   "<1%"           : TF_max > 0    & < 0.01
#'   "1%-10%"        : TF_max >= 0.01 & < 0.10
#'   ">=10%"         : TF_max >= 0.10
#'
#' @param vardict_filtered  Data frame — VarDict variants after standard somatic
#'                          filters (PoN, noise, VAF, alt-count, v_classes).
#'                          Typically mylist[['VarDict Variants Filt. 0.0050']].
#' @param mutect2_pass      Data frame — Mutect2 PASS variants after standard
#'                          somatic filters.
#'                          Typically mylist[['Mutect2 Variants PASS']].
#' @param chip_mchip        Data frame — classify_chip() output for M-CHIP.
#' @param chip_lpath        Data frame — classify_chip() output for L-CHIP
#'                          pathogenic.
#' @param chip_lput         Data frame — classify_chip() output for L-CHIP
#'                          putative.
#' @param sample_name       Character scalar — sample identifier, stored in
#'                          the returned tibble for downstream joining.
#'
#' @return A one-row tibble with columns:
#'   Sample                : sample identifier
#'   TF_max                : max VAF after CHIP and germline exclusion (NA if
#'                           undetectable); primary TF surrogate
#'   TF_p80                : 80th percentile VAF
#'   TF_mean               : arithmetic mean VAF
#'   TF_sd                 : standard deviation of VAFs
#'   TF_median             : median VAF
#'   TF_tier               : character tier label (see above)
#'   TF_n_somatic          : integer number of somatic variants in the pool
#'   TF_max_gene           : Hugo_Symbol of the max-VAF variant
#'   TF_max_hgvsp          : HGVSp_Short of the max-VAF variant
#'   TF_mutect2_concordant : logical — max-VAF variant present in Mutect2 PASS?
#'
#' @examples
#' tf_row <- estimate_tumor_fraction(
#'     vardict_filtered = mylist[['VarDict Variants Filt. 0.0050']],
#'     mutect2_pass     = mylist[['Mutect2 Variants PASS']],
#'     chip_mchip       = mylist[['VarDict CHIP M-CHIP']],
#'     chip_lpath       = mylist[['VarDict CHIP L-CHIP pathogenic']],
#'     chip_lput        = mylist[['VarDict CHIP L-CHIP putative']],
#'     sample_name      = sample
#' )
#' mylist[['Tumor Fraction']] <- tf_row

estimate_tumor_fraction <- function(vardict_filtered,
                                    mutect2_pass,
                                    chip_mchip,
                                    chip_lpath,
                                    chip_lput,
                                    sample_name) {
    
    # --- Build a position index for CHIP variants ----------------------------
    # Union of all three CHIP tables; use the same IDX key as the rest of the
    # pipeline (Chr_Pos_Ref_Alt) so exclusion is exact.
    .make_idx <- function(df) {
        if (nrow(df) == 0) return(character(0))
        paste(df$Chromosome, df$Start_Position,
              df$Tumor_Seq_Allele1, df$Tumor_Seq_Allele2, sep = "_")
    }
    
    chip_idx <- unique(c(
        .make_idx(chip_mchip),
        .make_idx(chip_lpath),
        .make_idx(chip_lput)
    ))
    
    # --- Exclude CHIP variants from the somatic pool -------------------------
    somatic <- vardict_filtered |>
        dplyr::mutate(
            .idx = paste(Chromosome, Start_Position,
                         Tumor_Seq_Allele1, Tumor_Seq_Allele2, sep = "_"),
            t_AF = suppressWarnings(as.numeric(t_AF))
        ) |>
        dplyr::filter(!.idx %in% chip_idx, !is.na(t_AF))
    
    # --- Exclude putative germline heterozygotes -----------------------------
    # Variants at VAF >= 0.35 are treated as constitutional (germline background)
    # following the same threshold used in 05_PanelOfNormals.R. This prevents
    # private germline heterozygotes that are rare in gnomAD from being picked
    # up as the max-VAF somatic variant and inflating the TF estimate.
    somatic <- somatic |>
        dplyr::filter(t_AF < 0.35)
    
    # --- Undetectable case ---------------------------------------------------
    if (nrow(somatic) == 0) {
        return(tibble::tibble(
            Sample                = sample_name,
            TF_max                = NA_real_,
            TF_p80                = NA_real_,
            TF_mean               = NA_real_,
            TF_sd                 = NA_real_,
            TF_median             = NA_real_,
            TF_tier               = "undetectable",
            TF_n_somatic          = 0L,
            TF_max_gene           = NA_character_,
            TF_max_hgvsp          = NA_character_,
            TF_mutect2_concordant = NA
        ))
    }
    
    # --- Descriptive statistics of the somatic VAF distribution -------------
    vafs     <- somatic$t_AF
    max_row  <- somatic |>
        dplyr::slice_max(order_by = t_AF, n = 1, with_ties = FALSE)
    
    tf_max    <- max_row$t_AF
    tf_p80    <- quantile(vafs, 0.80)
    tf_mean   <- mean(vafs)
    tf_sd     <- sd(vafs)
    tf_median <- median(vafs)
    
    # TF tier is driven by the max, following Husain et al. (2022).
    # The additional statistics characterise the VAF distribution:
    # a large gap between TF_max and TF_p80/TF_mean suggests the max-VAF
    # variant may be an outlier rather than representative of the bulk clone.
    tf_tier <- dplyr::case_when(
        tf_max >= 0.10 ~ ">=10%",
        tf_max >= 0.01 ~ "1%-10%",
        tf_max >  0    ~ "<1%",
        .default       = "undetectable"
    )
    
    # --- Mutect2 concordance check -------------------------------------------
    # Check whether the max-VAF VarDict variant is independently called as
    # PASS by Mutect2 at the same locus.
    max_idx <- paste(max_row$Chromosome, max_row$Start_Position,
                     max_row$Tumor_Seq_Allele1, max_row$Tumor_Seq_Allele2,
                     sep = "_")
    
    m2_concordant <- FALSE
    if (nrow(mutect2_pass) > 0) {
        m2_idx <- paste(mutect2_pass$Chromosome, mutect2_pass$Start_Position,
                        mutect2_pass$Tumor_Seq_Allele1,
                        mutect2_pass$Tumor_Seq_Allele2, sep = "_")
        m2_concordant <- max_idx %in% m2_idx
    }
    
    # --- Return one-row tibble -----------------------------------------------
    tibble::tibble(
        Sample                = sample_name,
        TF_max                = round(tf_max,    4),
        TF_p80                = round(tf_p80,    4),
        TF_mean               = round(tf_mean,   4),
        TF_sd                 = round(tf_sd,     4),
        TF_median             = round(tf_median, 4),
        TF_tier               = tf_tier,
        TF_n_somatic          = nrow(somatic),
        TF_max_gene           = max_row$Hugo_Symbol,
        TF_max_hgvsp          = max_row$HGVSp_Short,
        TF_mutect2_concordant = m2_concordant
    )
}

# Estimate Tumor Fraction via Kernel Density Estimation ------------------

#' Estimate tumor fraction from the mode of the somatic VAF distribution
#' using kernel density estimation (KDE)
#'
#' Non-parametric alternative to the max-VAF surrogate. The VAF distribution
#' of somatic variants above the noise floor is smoothed with a KDE and the
#' mode (highest-density peak) is identified. TF is estimated as 2 * mode,
#' assuming heterozygosity of the dominant somatic clone.
#'
#' This approach is complementary to estimate_tumor_fraction(): where the
#' max-VAF approach is sensitive to single outlier variants (flagged by
#' TF_max_outlier), the KDE mode reflects the bulk of the VAF distribution
#' and is more robust when enough variants are present. The two estimates
#' should broadly agree for well-behaved samples; a large discrepancy
#' suggests either a single dominant outlier (high TF_max, low TF_kde) or
#' a genuinely heterogeneous VAF distribution worth inspecting.
#'
#' The noise floor is estimated adaptively per sample as the
#' noise_floor_quantile of the VAF distribution (default 40th percentile).
#' This consistently separates the noise cluster from signal in low-TF
#' liquid biopsy samples without requiring a fixed threshold.
#'
#' The Sheather-Jones bandwidth selector (bw = "SJ") is used as it is
#' data-driven and generally performs well for unimodal and moderately
#' skewed distributions. For very sparse above-floor variant sets (<5
#' variants) the estimate is suppressed and NA returned.
#'
#' NOTE: This function is experimental. It is most reliable when
#' TF_n_somatic >= 15 and TF_max_outlier is TRUE (i.e. the max-VAF approach
#' is flagged as potentially unreliable). For typical low-TF samples with
#' few above-floor variants the KDE mode will sit just above the noise
#' floor and should not be over-interpreted.
#'
#' @param somatic               Data frame — same somatic pool used by
#'                              estimate_tumor_fraction(): CHIP-excluded,
#'                              germline-guarded (t_AF < 0.35), with numeric
#'                              t_AF column. Passed directly from the loop
#'                              after somatic pool construction.
#' @param noise_floor_quantile  Numeric in (0, 1) — quantile of the VAF
#'                              distribution used as the adaptive noise floor
#'                              cutoff (default 0.40).
#' @param min_above_floor       Integer — minimum number of variants above
#'                              the noise floor required to attempt KDE
#'                              (default 5). Returns NA silently otherwise.
#' @param min_variants          Integer — minimum total somatic variants
#'                              required before attempting KDE (default 15).
#'                              Returns NA silently otherwise.
#'
#' @return A one-row tibble with columns:
#'   TF_kde          : KDE mode-based TF estimate (2 * mode VAF); NA if
#'                     conditions not met or bandwidth selection fails
#'   TF_kde_floor    : adaptive noise floor used (the noise_floor_quantile
#'                     percentile of the VAF distribution)
#'   TF_kde_n_above  : number of variants above the noise floor used for
#'                     the KDE
#'
#' @examples
#' # In the per-sample loop, after somatic pool construction:
#' tf_kde_row <- estimate_tf_kde(somatic_loop)
#' tf_row     <- dplyr::bind_cols(tf_row, tf_kde_row)

estimate_tf_kde <- function(somatic,
                            noise_floor_quantile = 0.40,
                            min_above_floor      = 5L,
                            min_variants         = 15L) {
    
    na_result <- tibble::tibble(
        TF_kde         = NA_real_,
        TF_kde_floor   = NA_real_,
        TF_kde_n_above = NA_integer_
    )
    
    # --- Guard: enough total variants? --------------------------------------
    if (nrow(somatic) < min_variants) {
        message(sprintf(
            "estimate_tf_kde(): %d variants < min_variants (%d). Returning NA.",
            nrow(somatic), min_variants
        ))
        return(na_result)
    }
    
    vafs <- somatic$t_AF
    
    # --- Adaptive noise floor -----------------------------------------------
    # The noise_floor_quantile percentile of the full VAF distribution
    # separates the noise cluster from signal without requiring a fixed
    # threshold. In low-TF samples this consistently sits at ~0.006-0.009.
    noise_floor  <- quantile(vafs, noise_floor_quantile)
    vafs_signal  <- vafs[vafs > noise_floor]
    n_above      <- length(vafs_signal)
    
    # --- Guard: enough above-floor variants? --------------------------------
    if (n_above < min_above_floor) {
        message(sprintf(
            "estimate_tf_kde(): only %d variants above noise floor (%.4f). Returning NA.",
            n_above, noise_floor
        ))
        return(na_result)
    }
    
    # --- KDE with Sheather-Jones bandwidth ----------------------------------
    # bw = "SJ" is data-driven and avoids over-smoothing relative to the
    # default Silverman rule-of-thumb, which tends to over-smooth right-
    # skewed distributions like this one.
    dens <- tryCatch(
        density(vafs_signal,
                bw   = "SJ",
                from = noise_floor,
                to   = 0.35,
                n    = 512L),
        error = function(e) {
            message("estimate_tf_kde(): density() failed: ", conditionMessage(e))
            NULL
        }
    )
    
    if (is.null(dens)) return(na_result)
    
    # Mode of the density above the noise floor
    tf_half <- dens$x[which.max(dens$y)]
    
    tibble::tibble(
        TF_kde         = round(2 * tf_half, 4),
        TF_kde_floor   = round(noise_floor, 4),
        TF_kde_n_above = as.integer(n_above)
    )
}


