# Credentials--------------------------------------------------------------

#
# Creation of panel of normals for liquid biopsy workflow
# Author: Axel Künstner

# SPDX-License-Identifier: MIT
# Copyright (c) 2026 Axel Künstner
# Part of the LION panel analysis code: https://github.com/kunstner/UKSH_liquidbiopsy
# For research use only; not validated as a medical device.

# Libraries ---------------------------------------------------------------

pacman::p_load(tidyverse, openxlsx2, data.table)

# Panel genes
gene_list <- read_delim(file = '../data/Gene_list.txt', 
                        show_col_types = FALSE) |> 
    dplyr::mutate(Gene = case_when(
        Gene == 'FAM46C (TENT5C)' ~ 'TENT5C',
        .default = Gene
    )) |> 
    arrange(Gene)

# File with columns Batch and ID
PoNFile <- '../PoN_Samples.xlsx'
samples <- openxlsx2::read_xlsx(file = PoNFile)

list_pon <- list()

# Data --------------------------------------------------------------------

samples <- samples |> 
    dplyr::mutate(ID = gsub('cf', 'cn', ID)) |> 
    dplyr::mutate(Path = paste0('../', Batch, '/fgbio_results/', ID, '/vep.VarDict.MAF.gz')) 

dat.variants <- NULL
for( i in 1:nrow(samples) ) {
    # Transcript_Id not always present
    tmp <- read_delim(file = samples$Path[i], 
               skip = 1, show_col_types = FALSE, na = '#N/A',
               name_repair = 'unique_quiet') |> 
        dplyr::mutate(Tumor_Sample_Barcode = samples$ID[i]) |> 
        dplyr::mutate(t_AF = t_alt_count/t_depth) |> 
        dplyr::select( -ends_with('77'), -ends_with('136') ) |> 
        dplyr::filter(t_alt_count > 2) |> 
        dplyr::select(Tumor_Sample_Barcode,
                      Hugo_Symbol, Chromosome, Start_Position,
                      Ref = Reference_Allele, Alt = Tumor_Seq_Allele2,
                      HGVSc, HGVSp_Short, 
                      Variant_Type, Variant_Classification,
                      Existing_variation,
                      t_ref_count, t_alt_count, t_depth, t_AF,
                      MAX_AF, starts_with('gnomAD')) |> 
        dplyr::mutate(Background = case_when(
            t_AF >= 0.35 ~ TRUE,
            .default = FALSE
        )) |> 
        dplyr::mutate(Technical_background = case_when(
            t_AF < 0.35 ~ TRUE,
            .default = FALSE
        ))
    dat.variants <- rbindlist( list(dat.variants, tmp), fill = TRUE)
}

# Decode URL-encoded characters in HGVS notation
dat.variants <- dat.variants |>
    dplyr::mutate(dplyr::across(dplyr::any_of(c("HGVSc", "HGVSp", "HGVSp_Short")),
                                \(x) stringr::str_replace_all(x, "%3D", "=")))


list_pon[['Background']] <- dat.variants |> 
    dplyr::filter(Background == TRUE) |> 
    dplyr::select(-Tumor_Sample_Barcode, -t_ref_count, -t_alt_count, 
                  -t_depth, -t_AF, 
                  -Background, -Technical_background) |> 
    group_by( across( everything() ) ) |> 
    summarise(Pop_n = n(), .groups = 'drop') |> 
    dplyr::mutate(Pop_freq = round(Pop_n/nrow(samples), 4)) |> 
    dplyr::filter( Pop_n >= 2)
nrow(list_pon[['Background']]) 

# Implementation vaf noise ------------------------------------------------

# vaf noise
# (1) remove highest vaf per gene
# (2) sort vafs; max(vaf) * noise_factor -> threshold for noise vaf_noise
# (3) vaf_patient < vaf_noise -> FLAG as noise

noise_factor <- 1.25

dat.variants <- dat.variants |> 
    dplyr::mutate(IDX = paste(Hugo_Symbol, Chromosome, Start_Position, Ref, Alt, sep = "_")) |> 
    dplyr::mutate(IDX = factor(IDX))

df.noise <- data.frame(Hugo_Symbol = NULL, 
                       Chromosome = NULL, 
                       Start_Position = NULL, 
                       Ref = NULL, 
                       Alt = NULL, 
                       n_samples = NULL,
                       vaf_noise = NULL)

for( i in levels(dat.variants$IDX ) ) {
    tmp <- dat.variants |> 
        dplyr::filter(IDX == i) |> 
        dplyr::filter(t_alt_count > 2) |>
        dplyr::select(Tumor_Sample_Barcode, Hugo_Symbol, 
                      Chromosome, Start_Position, Ref, Alt, t_AF) |> 
        arrange(t_AF) |> # arrange data in ascending order
        dplyr::filter(row_number() <= n()-1) # remove last row -> highest value
    if( nrow(tmp >= 1) ) {
        vaf_noise <- min( 1, noise_factor * max(tmp$t_AF) )
        df.noise <- rbind(df.noise, data.frame(Hugo_Symbol = tmp$Hugo_Symbol[1], 
                   Chromosome = tmp$Chromosome[1], 
                   Start_Position = tmp$Start_Position[1], 
                   Ref = tmp$Ref[1], 
                   Alt = tmp$Alt[1], 
                   n_samples = nrow(tmp)+1,
                   vaf_noise = vaf_noise) )
    } else {
        # skip
    }
    rm(tmp)
}

list_pon[['Noise']] <- df.noise
nrow(list_pon[['Noise']]) 

# Export data -------------------------------------------------------------

# WriteXLS::WriteXLS(x = list_pon, 
#                    ExcelFileName = 'Panel_of_normal_fgbio.xlsx', 
#                    SheetNames = c('Background', 'Technical background', 'Noise'), 
#                    AdjWidth = T, BoldHeaderRow = T, FreezeRow = 1)

wb <- openxlsx2::wb_workbook()

for (sheet in names(list_pon)) {
    wb <- wb |>
        openxlsx2::wb_add_worksheet(sheet = sheet) |>
        openxlsx2::wb_add_data(sheet = sheet, x = list_pon[[sheet]]) |>
        openxlsx2::wb_freeze_pane(sheet = sheet, first_row = TRUE) |>
        openxlsx2::wb_add_font(
            sheet = sheet,
            dims  = openxlsx2::wb_dims(rows = 1, cols = 1:ncol(list_pon[[sheet]])),
            bold  = TRUE
        )
}

openxlsx2::wb_save(wb, file = 'Panel_of_normal_fgbio.xlsx', overwrite = TRUE)
