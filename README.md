# LION / HiLiB: R scripts for cfDNA variant summarization

R scripts for the downstream analysis of liquid biopsy (cfDNA) sequencing data of the **LION panel** (109-gene hybrid-capture panel with UMIs). They take the output of the fgbio-based preprocessing and variant calling workflow and produce annotated, filtered per-sample variant tables, CHIP classification, surrogate tumor fraction (TF) estimates, QC summaries and plots.

The upstream data preparation workflow (**`fgbio_v1_GRCh38`**: alignment, UMI consensus calling, VarDict/Mutect2 calling, VEP/vcf2maf annotation, fusion detection) is available at:
<https://github.com/kunstner/UKSH_liquidbiopsy>

The R scripts in this repository are the second stage of the analysis and expect the output of `fgbio_v1_GRCh38` (see [Upstream workflow](#upstream-workflow-fgbio_v1_grch38)).

> **Status:** accompanying code of the LION panel manuscript (Feierabend, Künstner et al.). Please cite the publication when using this code (see [Citation](#citation)).

---

## Overview

| Script | Purpose |
|---|---|
| `00_helper.R` | Helper functions: variant processing, PoN/noise annotation, drug annotation, filtering, QC and VAF plots, gene heatmap, tumor fraction estimation (max-VAF and experimental KDE) |
| `01_chip_classifier.R` | CHIP classification (M-CHIP, L-CHIP pathogenic, L-CHIP putative) including a binomial test to flag putative germline variants |
| `05_PanelOfNormals.R` | Builds the Panel of Normals (PoN) from healthy donor samples (recurrent background variants and per-variant VAF noise thresholds) |
| `10_summarize_fgbio.R` | Main per-batch script: reads VarDict/Mutect2 MAFs and QC files, applies filters, classifies CHIP, estimates TF, writes Excel tables, QC summaries and plots |

**Run order:** `05_PanelOfNormals.R` (once, or whenever the PoN samples change) -> `10_summarize_fgbio.R` (per batch). `00_helper.R` and `01_chip_classifier.R` are sourced by `10_summarize_fgbio.R`.

---

## Upstream workflow (`fgbio_v1_GRCh38`)

Raw sequencing data are processed per sample by the bash workflow `fgbio_v1_GRCh38` (repository: <https://github.com/kunstner/UKSH_liquidbiopsy>), designed for libraries prepared with the SureSelect XT HS2 kit (read structure `3M2S+T 3M2S+T`) and the human reference genome **GRCh38**.

Main steps:

1. UMI extraction (fgbio `ExtractUmisFromBam`), quality/adapter filtering (fastp)
2. Alignment (bwa-mem2) and merging with the unmapped BAM to retain UMI tags (Picard)
3. Duplicate marking (Picard `MarkDuplicates`) for internal QC
4. UMI-based error correction: read family grouping (`GroupReadsByUmi`, paired strategy), consensus calling (`CallMolecularConsensusReads`), realignment and `FilterConsensusReads`
5. Clipping of overlapping read pairs (`ClipBam`)
6. QC metrics (Picard `CollectTargetedPcrMetrics`, insert sizes, duplex/family size metrics)
7. Variant calling with **VarDict** (minimum VAF set via `AF_THR`) and **Mutect2**, left-alignment, VEP annotation and conversion to MAF (vcf2maf)
8. Fusion detection with FACTERA (structural variant calling with svict is present but disabled)

Per sample, the workflow exports the files listed under [Input data](#from-the-upstream-workflow) (`vep.VarDict.MAF.gz`, `vep.mutect2.MAF.gz`, `*_pcr_metrics.txt.gz`, `*_stats.insert.txt.gz`, `fusions_factera.txt.gz`, `duplex_metrics.*`, `family_size_hist.txt.gz`).

Before running the workflow, set the following variables (the script contains lab-specific defaults for a SLURM/HPC environment that need to be adapted):

| Variable | Meaning |
|---|---|
| `TARGET` | Base directory containing `001_resources`, `002_fastq` and `003_results` |
| `BATCH`, `SAMPLEID` | Batch and sample identifiers |
| `FASTQR1a/b`, `FASTQR2a/b` | Forward/reverse FASTQ files (two lanes are concatenated per sample) |
| `AF_THR` | Minimum VAF passed to VarDict (`0.0005` = 0.05%) |
| `MINFAMSIZE` | Minimum UMI family size for consensus calling and filtering (`2` in the manuscript) |
| `SCRATCH` | Temporary directory |

Required resources (paths defined in the script header): fgbio, fastp, bwa-mem2 with GRCh38 index, GATK, Picard, VarDict, VEP (Singularity image) with cache and plugin annotations (AlphaMissense, CADD, dbNSFP, gnomAD), vcf2maf, samtools, bcftools, FACTERA, the panel BED file, and the GRCh38 reference and GATK resource files (gnomAD germline resource, 1000G PoN, ExAC common variants).

> **Note:** The results of the workflow are written to `<TARGET>/003_results/<sample>/`. The R scripts expect them under `<BATCH>/fgbio_results/<sample>/`. Arrange or link the files accordingly.

---

## Requirements

- R (>= 4.1; the native pipe `|>` and `\(x)` lambda syntax are used)
- R packages (loaded via `pacman::p_load()`, which installs missing packages at runtime):
  `tidyverse`, `openxlsx2`, `patchwork`, `data.table`

For reproducible analyses we recommend pinning package versions, e.g. with `renv`.

---

## Input data

### From the upstream workflow

For each batch `<BATCH>` (e.g. `batch_a014`), the scripts expect the following under `<BATCH>/fgbio_results/<sample>/`:

| File | Content |
|---|---|
| `vep.VarDict.MAF.gz` | VarDict calls, VEP-annotated, MAF format |
| `*mutect2.MAF.gz` | Mutect2 calls, VEP-annotated, MAF format |
| `01_pcr_metrics.txt.gz` | QC metrics: mapping |
| `03_consensus_pcr_metrics.txt.gz` | QC metrics: after error correction |
| `04_pcr_metrics.txt.gz` | QC metrics: after remapping of clipped reads |
| `04_stats.insert.txt.gz` | Insert size statistics |
| `fusions_factera.txt.gz` | Fusion calls (FACTERA) |

The number of files per type must be identical across all samples of a batch; the script stops otherwise.

### Panel and reference files

| File | Content |
|---|---|
| `data/Gene_list.txt` | Panel gene list (tab-delimited; columns `Gene` and `Inclusion` are required) |
| `DB/` | CHIP gene lists (see below) |
| `Panel_of_normal_fgbio.xlsx` | Panel of Normals, created by `05_PanelOfNormals.R` (sheets `Background` and `Noise`) |
| `PoN_Samples.xlsx` | Sample sheet for the PoN (columns `ID`, `Batch`); required only for `05_PanelOfNormals.R`. **Not included: contains sample identifiers.** |

### CHIP gene lists (`DB/`)

The CHIP classification uses three gene lists (SARAH CH gene lists, v1.2):

- `M-CHIP_Genes_modified_v1.2.tsv`
- `L-CHIP_Genes_pathogenic_modified_v1.2.tsv`
- `L-CHIP_Genes_putative_modified_v1.2.tsv`

They are stored in the `DB/` folder. Make sure the constants `DBPATH_MCHIP`, `DBPATH_LCHIP_PATH` and `DBPATH_LCHIP_PUT` at the top of `10_summarize_fgbio.R` point to these files.

Position criteria in the `CHIP_Variants` columns are interpreted as **amino acid (protein) positions** by default and are matched against `HGVSp_Short`. Only tokens with an explicit `c.` prefix are treated as cDNA positions.

### DrugBank (`DRUGDB`): must be obtained separately

> **Important:** The DrugBank data are **not** part of this repository. DrugBank is distributed under a license that does not allow redistribution. To use the drug annotation, you must download the data yourself from <https://go.drugbank.com> (a license/account is required, depending on your use case) and place them in a local folder of your choice.

The scripts were developed with DrugBank version 5.1.12 and expect these two files:

- `all drug links.csv`
- `drugbank_all_target_polypeptide_ids/pharmacologically_active.csv`

Adjust the corresponding paths in `10_summarize_fgbio.R` (section *Data DrugBank*) to your local copy. Without these files the script will stop at that step.

---

## Usage

### 1. Build the Panel of Normals

Requires variant calls of healthy donor samples and a sample sheet (`PoN_Samples.xlsx`).

```bash
Rscript 05_PanelOfNormals.R
```

Output: `Panel_of_normal_fgbio.xlsx`

Variants are assigned to the PoN *Background* if they occur in at least two healthy donors with VAF >= 35% (constitutional background). Per-variant VAF noise thresholds are calculated as `1.25 x` the highest VAF observed in the healthy donors after removing the single highest value.

### 2. Summarize a batch

Interactive: set `BATCHID` in `10_summarize_fgbio.R` and source the script.

Command line (overrides `BATCHID`):

```bash
Rscript 10_summarize_fgbio.R batch_a014
```

Several batches:

```bash
for i in $(seq -w 1 14); do
    Rscript 10_summarize_fgbio.R "batch_a0${i}"
done
```

All paths in the scripts are **relative** to the working directory (`../<BATCH>/...`, `../data/...`). Start the scripts from the directory in which they are located, or adjust the paths in the setup section of `10_summarize_fgbio.R`.

---

## Output

Written to `<BATCH>/fgbio_tables/`:

| File | Content |
|---|---|
| `Results_<sample>_Table.xlsx` | Per-sample workbook: QC, tumor fraction, VarDict variants (unfiltered and filtered at 0.5%, 0.25% and 0.05% VAF), VarDict and Mutect2 CHIP classifications, Mutect2 variants (unfiltered and PASS), fusions |
| `Results_<sample>_heatmap.pdf` | Mutated genes per sample (combined, VarDict, Mutect2) |
| `Results_<sample>_vaf.pdf` | VAF distribution |
| `QC_Summary_<BATCH>.xlsx` | Batch-level QC (`QC_Summary`, `QC_Detailed`, `QC_All_Raw`) and tumor fraction summary (`TF_Summary`) |
| `QC_Summary_<BATCH>_plots.pdf` | Coverage, read counts and insert sizes |
| `QC_Advanced_<BATCH>_plots.pdf` | On-target efficiency, coverage uniformity, read retention, coverage variability |

---

## Method notes

- **Filtering:** population frequency (`MAX_AF` <= 0.001), minimum alt read count (> 2), variant classes of interest, exclusion of PoN background variants; the per-variant PoN noise threshold is annotated as `t_vaf_above_noise`.
- **CHIP classification:** gene- and variant-specific matching against the three CHIP lists; the binomial test (H0: VAF = 0.5 and H0: VAF = 1.0) flags variants consistent with germline heterozygosity or homozygosity (`CHIP_binom_flag`). Known high-VAF CHIP exceptions are always called somatic.
- **Tumor fraction:** surrogate TF following Husain et al. (2022, *JCO Precis Oncol*): the maximum VAF among somatic variants after exclusion of CHIP variants and putative germline variants (VAF >= 35%). Additional VAF distribution statistics (p80, mean, SD, median) and an experimental KDE-based estimate (`TF_kde`) are reported. Interpret TF estimates with caution in samples with few somatic variants.
- **VAF thresholds:** the thresholds used for the filtered tables are reporting thresholds of the bioinformatic workflow. They are not validated limits of detection (see the manuscript for the analytical validation).

---

## Notes and limitations

- The PoN is built from VarDict calls of healthy donors and applied to both callers.
- The scripts were developed for the LION panel and the fgbio-based workflow; use with other panels requires adjustments of the gene list, PoN and CHIP lists.
- The code is provided for research use and has not been validated as an in-vitro diagnostic.

---

## Citation

If you use this code, please cite:

> Feierabend S, Künstner A, et al. *A liquid biopsy-centered, pan-cancer, open next generation sequencing panel to support clinical decision-making (LION panel).* [journal, year, DOI: TODO]

Code archive: [Zenodo DOI: TODO]

The CHIP classification is adapted from the CHIP pipeline by L. Schawe; the ACMG classification function is taken from the MIRACUM pipeline. Please acknowledge these sources accordingly.

---

## License

[TODO: add license, e.g. MIT or GPL-3.0]

Third-party resources (DrugBank, CHIP gene lists) are subject to their own licenses and terms of use.

---

## Contact

Axel Künstner, University of Lübeck
