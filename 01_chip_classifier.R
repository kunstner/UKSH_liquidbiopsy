# Credentials--------------------------------------------------------------

#
# CHIP Classification Functions for cfDNA Analysis Pipeline
# Author: Axel Künstner
# Adapted from CHIP-Pipeline by @lschawe (2025-07-22)

# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Axel Künstner
# Part of the LION panel analysis code: https://github.com/kunstner/UKSH_liquidbiopsy
# For research use only; not validated as a medical device.

#
# Functions
#   load_chip_db          : load and pre-process a CHIP gene database file
#   .parse_cdna            : (internal) extract cDNA start/stop from CDS_position
#   .parse_aa_pos           : (internal) extract leading AA position from HGVSp_Short
#   .parse_exons            : (internal) extract exon numbers from Exon_Number field
#   .check_range            : (internal) check if a position falls within any given range
#   .check_mutation         : (internal) match a variant against per-type CHIP criteria;
#                            returns match-type string or NA_character_
#   .met_criteria            : (internal) dispatch variant to the correct criteria column;
#                            returns named list(match_type, criteria) or NA values
#   .chip_binom_test         : (internal) binomial test to flag putative germline variants
#   classify_chip            : main entry point — annotates variants with CHIP class
#
# Position convention (2026 rebuild):
#   Range and exact-residue criteria in the CHIP_Variants columns default to
#   AMINO ACID (protein) position, matched against HGVSp_Short. A token is
#   only treated as a cDNA/codon position when explicitly prefixed with
#   "c." (e.g. "c.741-791"). This was previously reversed (defaulted to
#   cDNA position), which silently suppressed matches for every gene using
#   bare or "p."-prefixed AA ranges (DNMT3A, TET2, IDH1/2, BRAF, CALR, CBL,
#   KIT, KMT2D, NOTCH1/2, NPM1, PHF6, PTPN11, RUNX1, SETBP1, SF3B1, SH2B3,
#   SRSF2, STAT3, GATA2/3, EZH2, FLT3, MPL, KRAS, NRAS, CREBBP).
#
# Output columns added by classify_chip():
#   CHIP_type          : "M-CHIP" | "L-CHIP_pathogenic" | "L-CHIP_putative"
#   CHIP_criteria      : human-readable matching rule from CHIP_Variants field of the db
#   CHIP_match_type    : computational evidence — "mutation_type_all" | "all" |
#                        "exon_match" | "cdna_range" | "aa_range" | "aa_match"
#   CHIP_binom_p       : p-value from exact binomial test (H0: VAF = 0.5)
#   CHIP_binom_CI_low  : lower bound of 95% CI on true VAF
#   CHIP_binom_CI_high : upper bound of 95% CI on true VAF
#   CHIP_binom_flag    : "somatic" | "germline_het" | "exception" | NA


# Libraries ---------------------------------------------------------------

pacman::p_load(tidyverse)


# Load CHIP Database ------------------------------------------------------

#' Load and pre-process a CHIP gene database TSV file
#'
#' @param db_path   Path to the TSV file (passed as a named constant, e.g. DBPATH_MCHIP)
#' @param chip_type One of "M-CHIP", "L-CHIP_pathogenic", "L-CHIP_putative"
#'
#' @return A tibble with columns: Gene, Mutation_type, Frameshift, Nonsense,
#'         Splice-Site, Missense, Deletion, Insertion, Delins, Transcript_ID
#'         Rows where CHIP_Variants == NA are dropped (no qualifying variants defined).
#'
#'         `na = "NA"` is used on read (not readr's default `c("", "NA")`)
#'         so a literal "NA" cell (category not applicable to this gene) stays
#'         distinguishable from a blank cell (category applicable, no residue
#'         restriction -> matches every variant of that class).
#'
#' @examples
#' chip_db_mchip  <- load_chip_db(DBPATH_MCHIP,      "M-CHIP")
#' chip_db_lpath  <- load_chip_db(DBPATH_LCHIP_PATH,  "L-CHIP_pathogenic")
#' chip_db_lput   <- load_chip_db(DBPATH_LCHIP_PUT,   "L-CHIP_putative")

load_chip_db <- function(db_path, chip_type) {
    
    allowed <- c("M-CHIP", "L-CHIP_pathogenic", "L-CHIP_putative")
    if (!chip_type %in% allowed) {
        stop(paste0(
            "'chip_type' must be one of: ",
            paste(allowed, collapse = ", "),
            ". Got: '", chip_type, "'"
        ))
    }
    
    db <- read_delim(file = db_path, delim = "\t", show_col_types = FALSE,
                     na = "NA") |>
        # Drop genes with no qualifying variants defined in the database
        dplyr::filter(!is.na(CHIP_Variants)) |>
        dplyr::mutate(
            # Trim stray whitespace on gene symbols (e.g. "SPRED2 ") so
            # dplyr::filter(chip_db, Gene == Hugo_Symbol) never silently fails
            Gene          = stringr::str_squish(Gene),
            # Derive a version-stripped transcript ID for MANE transcript matching
            # e.g. NM_015338.6 -> NM_015338
            Transcript_ID = stringr::str_remove(MANE_NM_v1.4, "\\.\\d+$"),
            CHIP_type     = chip_type,
            # Trim whitespace in every criteria column without disturbing NA
            dplyr::across(
                dplyr::any_of(c("Mutation_type", "Frameshift", "Nonsense",
                                "Splice-Site", "Missense", "Silent",
                                "Deletion", "Insertion", "Delins")),
                \(x) dplyr::if_else(is.na(x), x, stringr::str_squish(x))
            )
        ) |>
        dplyr::select(
            Gene, CHIP_type, CHIP_Variants, Transcript_ID, Mutation_type,
            Frameshift, Nonsense, `Splice-Site`, Missense,
            Silent, Deletion, Insertion, Delins
        )
    
    message(
        "Loaded ", chip_type, ": ",
        nrow(db), " genes with qualifying variants from '",
        basename(db_path), "'"
    )
    
    return(db)
}


# Internal helpers --------------------------------------------------------

#' Extract cDNA start and stop positions from a CDS_position string
#' e.g. "123-456" -> list(start=123, stop=456); "789" -> list(start=789, stop=NA)
#' @noRd
.parse_cdna <- function(cds) {
    if (is.na(cds) || cds == "") return(list(start = NA_real_, stop = NA_real_))
    range <- stringr::str_extract(cds, "[0-9]+-[0-9]+")
    if (!is.na(range)) {
        lim <- as.numeric(stringr::str_split(range, "-")[[1]])
        return(list(start = lim[1], stop = lim[2]))
    }
    start <- as.numeric(stringr::str_extract(cds, "[0-9]+"))
    list(start = start, stop = NA_real_)
}


#' Extract the leading amino-acid position from HGVSp_Short
#' e.g. "p.T125=" -> 125; "p.R337C" -> 337; "p.Thr337Cys" -> 337
#' @noRd
.parse_aa_pos <- function(x) {
    if (is.na(x) || x == "") return(NA_real_)
    as.numeric(stringr::str_extract(x, "[0-9]+"))
}


#' Extract exon numbers from a VEP Exon_Number string
#' e.g. "11/23" -> "11"; "11-12/23" -> c("11", "12")
#' @noRd
.parse_exons <- function(x) {
    if (is.na(x) || x == "") return(NA_character_)
    # VEP format: exon_number/total, keep only the left side
    x <- stringr::str_extract(x, "^[^/]+")
    parts    <- unlist(stringr::str_split(x, "-"))
    ranges   <- stringr::str_extract_all(parts, "[0-9]+-[0-9]+")[[1]]
    nums     <- unlist(stringr::str_extract_all(parts, "[0-9]+"))
    if (length(ranges) > 0) {
        for (r in ranges) {
            lim  <- as.numeric(stringr::str_split(r, "-")[[1]])
            nums <- unique(c(nums, as.character(seq(lim[1], lim[2]))))
        }
    }
    as.character(unique(nums))
}


#' Check whether a position (start, stop) falls within any of the
#' pipe-delimited numeric ranges in the criteria string
#' @noRd
.check_range <- function(start, stop, ranges) {
    if (is.na(start)) return(FALSE)
    if (is.na(stop)) stop <- start   # point variant: treat as the single-position interval [start, start]
    for (r in ranges) {
        vals <- as.numeric(stringr::str_split(r, "-")[[1]])
        if (length(vals) == 2 &&
            start <= vals[2] &&
            stop  >= vals[1]) return(TRUE)
    }
    FALSE
}


#' Core matching logic for a single mutation-type column
#'
#' Evaluates whether a variant matches the criteria stored in one column
#' (e.g. Frameshift, Missense) of the chip_db row.
#'
#' Position convention: AMINO ACID (protein) position is the default for
#' every range and exact-residue criterion, matched against HGVSp_Short.
#' A token is treated as cDNA/codon position ONLY when explicitly prefixed
#' with "c." (e.g. "c.741-791"); bare numbers ("290-912") and "p."-prefixed
#' numbers ("p.618-836") are always protein position.
#'
#' Matching hierarchy:
#'   1. NA criteria             -> not applicable to this gene: NA
#'   2. blank criteria          -> applicable, no restriction: "all"
#'   3. "all" token              -> unconditional pass: "all"
#'   4. transcript guard         -> position-specific criteria on a non-MANE
#'                                 transcript: NA
#'   5. "exon" keywords          -> match parsed exon number: "exon_match"
#'   6. "c." ranges               -> cDNA/codon range overlap: "cdna_range"
#'   7. AA ranges (default)      -> protein position range overlap: "aa_range"
#'   8. exact AA tokens           -> boundary-anchored substring match in
#'                                 HGVSp_Short: "aa_match"
#'
#' @return Character string ("all"|"exon_match"|"cdna_range"|"aa_range"|
#'         "aa_match") or NA_character_ if no match.
#' @noRd
.check_mutation <- function(criteria, aa_short, exon,
                            cdna_start, cdna_stop,
                            mane, transcript_id) {
    
    # 1. Literal NA on file -> category not applicable to this gene
    if (is.na(criteria)) return(NA_character_)
    
    # 2. Applicable category, no residue restriction on file (blank/whitespace
    #    cell) -> every variant of this mutation class counts
    if (stringr::str_squish(criteria) == "") return("all")
    
    tokens <- unlist(stringr::str_split(criteria, "\\|"))
    tokens <- stringr::str_squish(tokens)
    tokens <- tokens[tokens != ""]
    if (length(tokens) == 0) return("all")
    
    # Expand compressed multi-substitution shorthand, e.g. "R2832H|C" meaning
    # "R2832H or R2832C": a bare orphaned single-letter residue token (no
    # position of its own) inherits the original-residue+position prefix
    # from the most recent preceding token that carried one, so "C" becomes
    # "R2832C" before any matching is attempted. Without this, such tokens
    # are effectively unreachable (a bare letter almost never appears next
    # to non-digit characters in real 1-letter HGVSp notation) and carry a
    # small spurious-match risk if 3-letter codes ever appear.
    current_prefix <- NA_character_
    for (i in seq_along(tokens)) {
        tok <- tokens[i]
        if (stringr::str_detect(tok, stringr::regex("^all$", ignore_case = TRUE))) {
            next
        }
        if (stringr::str_detect(tok, "^[A-Z]$")) {
            if (!is.na(current_prefix)) {
                tokens[i] <- paste0(current_prefix, tok)
            }
        } else {
            prefix_match <- stringr::str_match(tok, "^([A-Za-z]+[0-9]+)")
            if (!is.na(prefix_match[1, 2])) {
                current_prefix <- prefix_match[1, 2]
            }
        }
    }
    
    if (!is.na(aa_short)) aa_short <- stringr::str_squish(aa_short)
    
    # 3. "all" token -> unconditional pass
    if (any(stringr::str_detect(tokens, stringr::regex("^all$", ignore_case = TRUE)))) {
        return("all")
    }
    
    # Transcript guard for position-specific criteria: if we cannot confirm
    # this variant is annotated on the MANE transcript, position/range
    # criteria cannot be safely evaluated (numbering may differ under a
    # non-MANE transcript) -- fail the guard rather than pass through.
    if (is.na(mane) || is.na(transcript_id) || mane != transcript_id) {
        return(NA_character_)
    }
    
    # 5. Exon matching
    if (any(stringr::str_detect(tokens, stringr::regex("exon", ignore_case = TRUE)))) {
        crit_nums <- unlist(stringr::str_extract_all(tokens, "[0-9]+(?:-[0-9]+)?"))
        exon_nums <- .parse_exons(exon)
        if (any(exon_nums %in% unlist(stringr::str_extract_all(crit_nums, "[0-9]+")))) {
            return("exon_match")
        }
    }
    
    # 6. Explicit cDNA/codon range -- ONLY tokens carrying a "c." prefix
    c_tokens <- tokens[stringr::str_detect(tokens, stringr::regex("^c\\.", ignore_case = TRUE))]
    if (length(c_tokens) > 0) {
        cdna_ranges <- unlist(stringr::str_extract_all(c_tokens, "[0-9]+-[0-9]+"))
        if (length(cdna_ranges) > 0 && .check_range(cdna_start, cdna_stop, cdna_ranges)) {
            return("cdna_range")
        }
    }
    
    # 7. Amino-acid position range -- DEFAULT for every other ranged token,
    #    bare ("618-836") or "p."-prefixed ("p.618-836")
    aa_tokens <- tokens[!stringr::str_detect(tokens, stringr::regex("^c\\.", ignore_case = TRUE))]
    aa_ranges <- unlist(stringr::str_extract_all(aa_tokens, "[0-9]+-[0-9]+"))
    if (length(aa_ranges) > 0) {
        aa_pos <- .parse_aa_pos(aa_short)
        if (!is.na(aa_pos)) {
            hit <- vapply(aa_ranges, function(r) {
                lim <- as.numeric(stringr::str_split(r, "-")[[1]])
                aa_pos >= lim[1] && aa_pos <= lim[2]
            }, logical(1))
            if (any(hit)) return("aa_range")
        }
    }
    
    # 8. Exact AA / HGVSp token matching (e.g. "R337", "T615"), boundary-
    #    anchored on digits so "R337" cannot match "R3370"
    if (!is.na(aa_short) && aa_short != "") {
        exact_tokens <- tokens[!stringr::str_detect(tokens, "-[0-9]")]
        for (tok in exact_tokens) {
            # Bare frameshift hotspot tokens (e.g. "P95fs") never carry the
            # replacement residue / stop distance that real HGVS frameshift
            # notation always includes (e.g. "p.P95Rfs*3"). For these,
            # match on original residue + position, followed by ANY short
            # residue code, then "fs" -- rather than a literal substring.
            fs_hotspot <- stringr::str_match(tok, "^([A-Za-z])([0-9]+)fs\\*?[0-9]*$")
            if (!is.na(fs_hotspot[1, 1])) {
                orig_aa <- fs_hotspot[1, 2]
                pos     <- fs_hotspot[1, 3]
                pat <- paste0("(?<![0-9A-Za-z])", stringr::str_escape(orig_aa), pos,
                              "[A-Za-z]{0,3}fs")
            } else {
                pat <- paste0("(?<![0-9])", stringr::str_escape(tok), "(?![0-9])")
            }
            if (stringr::str_detect(aa_short, pat)) return("aa_match")
        }
    }
    
    NA_character_
}


#' Dispatch a variant to the correct criteria column based on Variant_Classification
#' and evaluate it against chip_db
#'
#' @param gene           Hugo_Symbol of the variant
#' @param var_class      Variant_Classification (VEP/MAF convention)
#' @param aa_short       HGVSp_Short (already normalized via .normalize_hgvsp())
#' @param exon           Exon_Number (VEP format, e.g. "11/23")
#' @param cdna_start     Numeric cDNA start (from .parse_cdna)
#' @param cdna_stop      Numeric cDNA stop  (from .parse_cdna)
#' @param mane           Version-stripped transcript ID of the variant
#' @param chip_db        Pre-loaded chip_db tibble (from load_chip_db)
#'
#' @return Named list:
#'   $match_type : "mutation_type_all" | "all" | "exon_match" | "cdna_range" |
#'                 "aa_range" | "aa_match" | NA_character_
#'   $criteria   : CHIP_Variants string from the database row, or NA_character_
#' @noRd
.met_criteria <- function(gene, var_class, aa_short, exon,
                          cdna_start, cdna_stop, mane, chip_db) {
    
    no_match <- list(match_type = NA_character_, criteria = NA_character_)
    
    row <- dplyr::filter(chip_db, Gene == gene)
    if (nrow(row) == 0) return(no_match)
    row <- row[1, ]
    
    transcript_id <- row$Transcript_ID
    chip_criteria <- row$CHIP_Variants   # human-readable rule from the database
    
    # If Mutation_type contains "all", every functional variant qualifies
    mt_types <- unlist(stringr::str_split(row$Mutation_type, "\\|"))
    if (any(stringr::str_detect(mt_types, stringr::regex("^all$", ignore_case = TRUE)))) {
        return(list(match_type = "mutation_type_all", criteria = chip_criteria))
    }
    
    # Dispatch by VEP Variant_Classification
    vc <- var_class
    
    match_col <- function(col) {
        .check_mutation(row[[col]], aa_short, exon,
                        cdna_start, cdna_stop, mane, transcript_id)
    }
    
    match_type <- dplyr::case_when(
        stringr::str_detect(vc, "Frame_Shift")       ~ match_col("Frameshift"),
        stringr::str_detect(vc, "Nonsense_Mutation")  ~ match_col("Nonsense"),
        stringr::str_detect(vc, "Splice_Site|Splice_Region") ~ match_col("Splice-Site"),
        # Delins carry "delins" in HGVSp_Short but are classified as Missense by VEP
        stringr::str_detect(vc, "Missense_Mutation") &
            !is.na(aa_short) &
            stringr::str_detect(aa_short, stringr::regex("delins", ignore_case = TRUE))
        ~ match_col("Delins"),
        stringr::str_detect(vc, "Missense_Mutation")  ~ match_col("Missense"),
        stringr::str_detect(vc, "In_Frame_Ins")       ~ match_col("Insertion"),
        stringr::str_detect(vc, "In_Frame_Del")       ~ match_col("Deletion"),
        stringr::str_detect(vc, "Silent")             ~ match_col("Silent"),
        .default = NA_character_
    )
    
    if (is.na(match_type)) return(no_match)
    list(match_type = match_type, criteria = chip_criteria)
}


#' Binomial test to flag variants consistent with germline heterozygosity
#'
#' Tests H0: p(alt) = 0.5.  A significant result (p < p_threshold) means the
#' VAF is inconsistent with a germline het variant, i.e. the variant is likely
#' somatic / true CHIP.
#'
#' Known high-VAF CHIP exceptions (Vlasschaert et al. 2023, Blood) are always
#' flagged as somatic regardless of the test result.
#'
#' @param variants     Data frame with columns t_alt_count, t_depth, HGVSp_Short
#' @param p_threshold  Significance threshold (default 0.0001)
#'
#' @return variants with four new columns:
#'   CHIP_binom_p       : raw p-value from binom.test
#'   CHIP_binom_CI_low  : lower bound of 95% CI on the true VAF
#'   CHIP_binom_CI_high : upper bound of 95% CI on the true VAF
#'   CHIP_binom_flag    : "somatic" | "germline_het" | "exception" | NA
#' @noRd
.chip_binom_test <- function(variants, p_threshold = 0.0001) {
    
    # Known high-VAF CHIP exceptions (Vlasschaert et al. 2023, Blood 141:2214)
    highvaf_exceptions <- tibble::tibble(
        Hugo_Symbol = "TET2",
        HGVSp_Short = c("p.H1904R", "p.I1873T", "p.T1884A")
    )
    
    # Run both tests in one pass per variant
    binom_results <- mapply(
        function(ad, dp) {
            if (is.na(ad) || is.na(dp) || dp == 0) {
                return(list(
                    p_het  = NA_real_, p_hom  = NA_real_,
                    ci_low = NA_real_, ci_high = NA_real_
                ))
            }
            bt_het <- stats::binom.test(x = ad, n = dp, p = 0.5, alternative = "two.sided")
            bt_hom <- stats::binom.test(x = ad, n = dp, p = 1.0, alternative = "two.sided")
            list(
                p_het   = bt_het$p.value,
                p_hom   = bt_hom$p.value,
                ci_low  = bt_het$conf.int[1],
                ci_high = bt_het$conf.int[2]
            )
        },
        variants$t_alt_count, variants$t_depth,
        SIMPLIFY = FALSE
    )
    
    variants |>
        dplyr::mutate(
            CHIP_binom_p_het   = vapply(binom_results, `[[`, numeric(1), "p_het"),
            CHIP_binom_p_hom   = vapply(binom_results, `[[`, numeric(1), "p_hom"),
            CHIP_binom_CI_low  = round(vapply(binom_results, `[[`, numeric(1), "ci_low"),  4),
            CHIP_binom_CI_high = round(vapply(binom_results, `[[`, numeric(1), "ci_high"), 4),
            CHIP_binom_flag    = dplyr::case_when(
                # Known exceptions are always called somatic
                paste0(Hugo_Symbol, "_", HGVSp_Short) %in%
                    paste0(highvaf_exceptions$Hugo_Symbol, "_",
                           highvaf_exceptions$HGVSp_Short)  ~ "exception",
                is.na(CHIP_binom_p_het)                     ~ NA_character_,
                # Consistent with germline homozygous -> flag
                CHIP_binom_p_hom >= p_threshold             ~ "germline_hom",
                # Consistent with germline heterozygous -> flag
                CHIP_binom_p_het >= p_threshold             ~ "germline_het",
                # Inconsistent with both -> somatic CHIP
                .default                                    = "somatic"
            )
        )
}

# Main entry point --------------------------------------------------------

#' Classify variants as CHIP using a pre-loaded CHIP gene database
#'
#' This function is the direct replacement for:
#'   apply_variant_filters(..., ch_genes = gene_ch)
#'
#' It applies the full CHIP matching logic (gene membership, mutation type,
#' and positional/AA specificity) and annotates each passing variant with
#' its CHIP class, match evidence, and a binomial germline test result.
#'
#' Variants that do not meet CHIP criteria are dropped.
#'
#' @param variants            Data frame of variants from process_variants().
#'                            Must contain: Hugo_Symbol, Variant_Classification,
#'                            HGVSp_Short, Exon_Number, CDS_position,
#'                            Transcript_ID, t_alt_count, t_depth, t_AF.
#' @param chip_db             Pre-loaded CHIP database from load_chip_db().
#' @param vaf_threshold       Minimum VAF (default 0.005); set NULL to skip.
#' @param alt_count_threshold Minimum alt read count (default 2); set NULL to skip.
#' @param run_binom_test      Logical; run binomial germline test (default TRUE).
#' @param p_threshold         Significance threshold for binom test (default 0.0001).
#'
#' @return A filtered and annotated data frame with additional columns:
#'   CHIP_type          : "M-CHIP" | "L-CHIP_pathogenic" | "L-CHIP_putative"
#'   CHIP_criteria      : human-readable matching rule from CHIP_Variants in the db
#'   CHIP_match_type    : "mutation_type_all" | "all" | "exon_match" |
#'                        "cdna_range" | "aa_range" | "aa_match"
#'   CHIP_binom_p       : p-value from exact binomial test (H0: VAF = 0.5)
#'   CHIP_binom_CI_low  : lower bound of 95% CI on true VAF
#'   CHIP_binom_CI_high : upper bound of 95% CI on true VAF
#'   CHIP_binom_flag    : "somatic" | "germline_het" | "exception" | NA
#'
#' @examples
#' # In 10_summarize_fgbio.R — load once at top of script:
#' chip_db_mchip  <- load_chip_db(DBPATH_MCHIP,      "M-CHIP")
#' chip_db_lpath  <- load_chip_db(DBPATH_LCHIP_PATH,  "L-CHIP_pathogenic")
#' chip_db_lput   <- load_chip_db(DBPATH_LCHIP_PUT,   "L-CHIP_putative")
#'
#' # Inside the per-sample loop:
#' mylist[['VarDict CHIP M-CHIP']]            <- classify_chip(vardict_variants, chip_db_mchip)
#' mylist[['VarDict CHIP L-CHIP pathogenic']] <- classify_chip(vardict_variants, chip_db_lpath)
#' mylist[['Mutect2 CHIP M-CHIP']]            <- classify_chip(mutect2_variants, chip_db_mchip)

classify_chip <- function(variants,
                          chip_db,
                          vaf_threshold       = 0.005,
                          alt_count_threshold = 2,
                          run_binom_test      = TRUE,
                          p_threshold         = 0.0001) {
    
    # --- Input guards --------------------------------------------------------
    required_cols <- c(
        "Hugo_Symbol", "Variant_Classification", "HGVSp_Short",
        "Exon_Number", "CDS_position", "MANE_SELECT",
        "t_alt_count", "t_depth", "t_AF"
    )
    missing_cols <- setdiff(required_cols, colnames(variants))
    if (length(missing_cols) > 0) {
        stop(paste0(
            "classify_chip(): missing required columns: ",
            paste(missing_cols, collapse = ", ")
        ))
    }
    
    # --- Normalize gene symbols and HGVSp notation before any matching -------
    # Trims stray whitespace on Hugo_Symbol (mirrors Gene trimming in
    # load_chip_db()) and applies .normalize_hgvsp() defensively -- even if
    # process_variants() already normalized upstream, this keeps
    # classify_chip() correct when called standalone.
    out <- variants |>
        dplyr::mutate(
            Hugo_Symbol = stringr::str_squish(Hugo_Symbol),
            HGVSp_Short = .normalize_hgvsp(HGVSp_Short) |> stringr::str_squish()
        )
    
    # --- Pre-filter: restrict to genes present in this chip_db ---------------
    chip_genes <- unique(chip_db$Gene)
    
    out <- out |>
        dplyr::filter(Hugo_Symbol %in% chip_genes)
    
    if (nrow(out) == 0) {
        message("classify_chip(): no variants found for genes in this chip_db.")
        return(out)
    }
    
    # --- Ensure t_AF is numeric ----------------------------------------------
    # Mutect2 MAFs often carry '.' or empty strings in t_AF; fall back to
    # computing VAF from counts when the column is not cleanly numeric.
    out <- out |>
        dplyr::mutate(
            t_AF = suppressWarnings(as.numeric(t_AF)),
            t_AF = dplyr::case_when(
                !is.na(t_AF)  ~ t_AF,
                t_depth > 0   ~ t_alt_count / t_depth,
                .default      = NA_real_
            )
        )
    
    # --- Optional VAF / alt count thresholds ---------------------------------
    if (!is.null(vaf_threshold)) {
        out <- dplyr::filter(out, t_AF >= vaf_threshold)
    }
    if (!is.null(alt_count_threshold)) {
        out <- dplyr::filter(out, t_alt_count > alt_count_threshold)
    }
    
    # --- Strip transcript version for matching (e.g. NM_015338.6 -> NM_015338)
    # MANE_SELECT carries VEP's RefSeq MANE annotation, matching the
    # convention used in chip_db$Transcript_ID (derived from MANE_NM_v1.4).
    # Transcript_ID itself is Ensembl (ENST...) and is NOT comparable to the
    # database's RefSeq IDs -- using it here silently failed every
    # position-specific CHIP match.
    out <- out |>
        dplyr::mutate(
            .transcript_stripped = stringr::str_remove(MANE_SELECT, "\\.\\d+$")
        )
    
    # --- Parse cDNA positions once (vectorised via mapply) -------------------
    cdna_parsed     <- mapply(.parse_cdna, out$CDS_position, SIMPLIFY = FALSE)
    out$.cdna_start <- vapply(cdna_parsed, `[[`, numeric(1), "start")
    out$.cdna_stop  <- vapply(cdna_parsed, `[[`, numeric(1), "stop")
    
    # --- Apply CHIP criteria row-wise ----------------------------------------
    # .met_criteria() now returns list(match_type, criteria) per variant
    criteria_results <- mapply(
        .met_criteria,
        gene       = out$Hugo_Symbol,
        var_class  = out$Variant_Classification,
        aa_short   = out$HGVSp_Short,
        exon       = out$Exon_Number,
        cdna_start = out$.cdna_start,
        cdna_stop  = out$.cdna_stop,
        mane       = out$.transcript_stripped,
        MoreArgs   = list(chip_db = chip_db),
        SIMPLIFY   = FALSE
    )
    
    out <- out |>
        dplyr::mutate(
            CHIP_match_type = vapply(criteria_results, `[[`, character(1), "match_type"),
            CHIP_criteria   = vapply(criteria_results, `[[`, character(1), "criteria")
        ) |>
        # Keep only variants that matched a CHIP rule
        dplyr::filter(!is.na(CHIP_match_type)) |>
        dplyr::mutate(CHIP_type = unique(chip_db$CHIP_type)) |>
        # Drop internal working columns
        dplyr::select(-dplyr::starts_with("."))
    
    if (nrow(out) == 0) {
        message("classify_chip(): no variants passed CHIP criteria.")
        return(out)
    }
    
    # --- Binomial germline test ----------------------------------------------
    if (run_binom_test) {
        out <- .chip_binom_test(out, p_threshold = p_threshold)
    }
    
    return(out)
}
