# LION: R scripts for cfDNA variant summarization

R scripts for the downstream analysis of liquid biopsy (cfDNA) sequencing data of the **LION panel** (109-gene hybrid-capture panel with UMIs). They take the output of the fgbio-based preprocessing and variant calling workflow and produce annotated, filtered per-sample variant tables, CHIP classification, surrogate tumor fraction (TF) estimates, QC summaries and plots.

The upstream data preparation workflow (**`fgbio_v1_GRCh38`**: alignment, UMI consensus calling, VarDict/Mutect2 calling, VEP/vcf2maf annotation, fusion detection) is available in this repository as well.

The R scripts in this repository are the second stage of the analysis and expect the output of `fgbio_v1_GRCh38` (see [Upstream workflow](#upstream-workflow-fgbio_v1_grch38)).

> **Status:** accompanying code of the LION panel manuscript (Feierabend, Künstner et al.). Please cite the publication when using this code (see [Citation](#citation)).

---

## Overview

| Script | Purpose |
|---|---|
| `fgbio_v1_GRCh38` | Upstream bash workflow (per sample): UMI extraction, alignment, UMI consensus calling, QC metrics, VarDict/Mutect2 calling, VEP annotation, MAF conversion, fusion detection |
| `00_helper.R` | Helper functions: variant processing, PoN/noise annotation, drug annotation, filtering, QC and VAF plots, gene heatmap, tumor fraction estimation (max-VAF and experimental KDE) |
| `01_chip_classifier.R` | CHIP classification (M-CHIP, L-CHIP pathogenic, L-CHIP putative) including a binomial test to flag putative germline variants |
| `05_PanelOfNormals.R` | Builds the Panel of Normals (PoN) from healthy donor samples (recurrent background variants and per-variant VAF noise thresholds) |
| `10_summarize_fgbio.R` | Main per-batch script: reads VarDict/Mutect2 MAFs and QC files, applies filters, classifies CHIP, estimates TF, writes Excel tables, QC summaries and plots |

**Run order:** `fgbio_v1_GRCh38` (per sample) -> `05_PanelOfNormals.R` (once, or whenever the PoN samples change) -> `10_summarize_fgbio.R` (per batch). `00_helper.R` and `01_chip_classifier.R` are sourced by `10_summarize_fgbio.R`.

---

## Upstream workflow (`fgbio_v1_GRCh38`)

Raw sequencing data are processed per sample by the bash workflow `fgbio_v1_GRCh38`, designed for libraries prepared with the SureSelect XT HS2 kit (read structure `3M2S+T 3M2S+T`) and the human reference genome **GRCh38**.

The workflow runs with `set -eo pipefail`: it stops at the first failing step (including failures inside pipes). Large intermediate files (FASTQ, unmapped/aligned BAMs, read-family BAMs) are deleted as soon as the next step has succeeded, so the disk footprint per sample stays small. Only `ValidateSamFile` is non-fatal (findings are reported as a warning).

Main steps:

1. **Read preparation:** FASTQ to unmapped BAM (Picard `FastqToSam`), UMI extraction (fgbio `ExtractUmisFromBam`, UMIs stored in the `RX` tag), conversion back to FASTQ and quality/adapter filtering (fastp)
2. **Alignment:** bwa-mem2, streamed directly into a queryname sort (Picard `SortSam`, no intermediate SAM on disk), followed by `MergeBamAlignment` with the unmapped BAM to retain the UMI tags (coordinate-sorted output)
3. **Mapping QC:** `samtools flagstat`, Picard `CollectTargetedPcrMetrics` (`01_pcr_metrics.txt`)
4. **UMI-based error correction:** read family grouping (`GroupReadsByUmi`, paired strategy, 0 UMI edits, MAPQ >= 20), UMI/duplex metrics (`CollectDuplexSeqMetrics`), molecular consensus calling (`CallMolecularConsensusReads`, minimum family size `MINFAMSIZE`), realignment of consensus reads (bwa-mem2 + `MergeBamAlignment`) and consensus filtering (`FilterConsensusReads`); QC after error correction (`03_consensus_pcr_metrics.txt`)
5. **Clipping** of overlapping read pairs (`ClipBam`) and BAM validation (`ValidateSamFile`). The resulting consensus BAM is the input for all variant callers
6. **Final QC:** Picard `CollectTargetedPcrMetrics` (`04_pcr_metrics.txt`) and GATK `CollectInsertSizeMetrics` (`04_stats.insert.txt`)
7. **VarDict** (minimum VAF set via `AF_THR`, strand bias test, only non-filtered calls kept), left-alignment (GATK `LeftAlignAndTrimVariants`), VEP annotation and conversion to MAF (vcf2maf)
8. **Mutect2** (tumor-only, with 1000G PoN and gnomAD germline resource), read orientation model, contamination estimation, `FilterMutectCalls`, left-alignment, VEP annotation and conversion to MAF. All Mutect2 calls are kept in the MAF together with their `FILTER` value; selection of `PASS` calls is done in `10_summarize_fgbio.R`
9. **Fusion detection** with FACTERA
10. Clean-up and copy of the results to `003_results/<SAMPLEID>/`

Structural variant calling is currently not performed (placeholder step only). Duplicate marking (Picard `MarkDuplicates`) was part of earlier versions of the workflow and has been removed, as its output was used neither for variant calling nor by the R scripts.

VEP is run with `--mane --hgvs --numbers --canonical`, and `MANE_SELECT` is retained in both MAFs (`vcf2maf --retain-ann`). `MANE_SELECT`, `HGVSp_Short`, `Exon_Number` and `CDS_position` are required by the CHIP classifier (`01_chip_classifier.R`); `MAX_AF` is required for the population frequency filter.

### Exported files

Per sample, the workflow copies the following files to `<RESOURCES>/003_results/<SAMPLEID>/`:

| File | Content |
|---|---|
| `vep.VarDict.MAF.gz` | VarDict calls, VEP-annotated, MAF format |
| `vep.mutect2.MAF.gz` | Mutect2 calls (all FILTER values), VEP-annotated, MAF format |
| `01_pcr_metrics.txt.gz` | Targeted PCR metrics after mapping |
| `03_consensus_pcr_metrics.txt.gz` | Targeted PCR metrics after UMI consensus calling and filtering |
| `04_pcr_metrics.txt.gz` | Targeted PCR metrics of the final (clipped) consensus BAM |
| `04_stats.insert.txt.gz` | Insert size metrics of the final consensus BAM |
| `fusions_factera.txt.gz` | Fusion calls (FACTERA) |
| `duplex_metrics.*` | fgbio UMI/duplex metrics |
| `family_size_hist.txt.gz` | UMI family size histogram |

Further files (fastp reports, `01_stats_mapping.txt.gz`, `04_stats.insert.pdf`, VCFs, final BAM) remain in the per-sample working directory `<TARGET>/<SAMPLEID>/`.

### Variables

The following variables are not set inside the workflow and must be provided by the calling (starter/SLURM) script or exported in the environment:

| Variable | Meaning |
|---|---|
| `SAMPLEID` | Sample identifier; also the name of the per-sample working and results directory |
| `TARGET` | Working directory; a subdirectory `<TARGET>/<SAMPLEID>/` is created for all intermediate and output files |
| `FASTQR1a`, `FASTQR1b`, `FASTQR2a`, `FASTQR2b` | Forward (R1) and reverse (R2) FASTQ files; two lanes are concatenated per sample (for a single lane, leave the `b` variables empty) |
| `AF_THR` | Minimum VAF passed to VarDict (`0.0005` = 0.05%) |
| `MINFAMSIZE` | Minimum UMI family size for consensus calling and filtering (`2` in the manuscript) |
| `SCRATCH` | Scratch directory; temporary files are written to `$SCRATCH/tmp` |

The following variables are set in the script header and must be adapted to the local installation:

| Variable | Meaning |
|---|---|
| `RESOURCES` | Base directory containing `001_resources` (tools, reference, panel and GATK resource files) and `003_results` (exported results) |
| `ANNOTATIONS` | Directory containing the VEP cache (`vep/vep_data`) and the plugin/custom annotation files (`vep/annotations_GRCh38`) |
| `BEDFILE` | Panel target regions; use `DB/LION_Panel_GRCh38.bed` (see [Panel target regions](#panel-target-regions-db)) |

Example call (outside SLURM):

```bash
export SAMPLEID=sample_001 TARGET=/path/to/work SCRATCH=/path/to/scratch \
       AF_THR=0.0005 MINFAMSIZE=2 \
       FASTQR1a=/path/to/L001_R1.fastq.gz FASTQR1b=/path/to/L002_R1.fastq.gz \
       FASTQR2a=/path/to/L001_R2.fastq.gz FASTQR2b=/path/to/L002_R2.fastq.gz
bash fgbio_v1_GRCh38
```

### Required software and resources

- Java (for fgbio, Picard, GATK) and Singularity/Apptainer (for the VEP image), Perl (vcf2maf, FACTERA, VarDict helper scripts) and R (VarDict `teststrandbias.R`)
- fgbio 2.4.0, fastp 0.24.0, bwa-mem2, Picard 2.27.5, GATK 4.5.0.0, VarDictJava 1.8.3, samtools, VEP 111 (Singularity image), vcf2maf 1.6.22, FACTERA (requires BLAST+ and `twoBitToFa` in the `PATH`)
- GRCh38 reference (`Homo_sapiens_assembly38.fasta` with `.dict` and `.fai`) and bwa-mem2 index
- Gene annotation for FACTERA (`genes_GRCh38.bed`) and `hg38.2bit`
- GATK resource files: gnomAD germline resource (`af-only-gnomad.hg38.vcf.gz`), 1000G PoN (`1000g_pon.hg38.vcf.gz`), ExAC common variants (`small_exac_common_3.hg38.vcf.gz`)
- VEP cache (GRCh38) and annotations: AlphaMissense, CADD 1.6, dbNSFP 4.4c, gnomAD genomes r3.0 and gnomAD exomes r2.1.1 (liftover to GRCh38)

> **Note:** The workflow writes its results to `<RESOURCES>/003_results/<SAMPLEID>/`. The R scripts expect them under `../<BATCH>/fgbio_results/<SAMPLEID>/` (relative to the R scripts). Copy or link the per-sample result folders accordingly, grouped by batch.

---

## Requirements (R)

- R (>= 4.1; the native pipe `|>` and `\(x)` lambda syntax are used)
- `pacman` (packages are loaded via `pacman::p_load()`, which installs missing packages at runtime)
- R packages: `tidyverse`, `openxlsx2`, `patchwork`, `data.table` (the latter only for `05_PanelOfNormals.R`)

For reproducible analyses we recommend pinning package versions, e.g. with `renv`.

---

## Input data

### From the upstream workflow

For each batch `<BATCH>` (e.g. `batch_a001`), `10_summarize_fgbio.R` searches `../<BATCH>/fgbio_results/` recursively for:

| File | Content |
|---|---|
| `vep.VarDict.MAF.gz` | VarDict calls, VEP-annotated, MAF format |
| `vep.mutect2.MAF.gz` | Mutect2 calls, VEP-annotated, MAF format |
| `01_pcr_metrics.txt.gz` | QC metrics: mapping |
| `03_consensus_pcr_metrics.txt.gz` | QC metrics: after error correction |
| `04_pcr_metrics.txt.gz` | QC metrics: final consensus BAM (clipped reads) |
| `04_stats.insert.txt.gz` | Insert size statistics |
| `fusions_factera.txt.gz` | Fusion calls (FACTERA) |

The sample name is taken from the name of the folder containing `vep.VarDict.MAF.gz` (a prefix `fullrun_` is removed). Every sample of a batch must provide all seven files; the script checks that the number of files is identical for each file type and stops otherwise.

### Panel and reference files

| File | Content |
|---|---|
| `DB/Gene_list.txt` | Panel gene list (tab-delimited; columns `Gene` and `Inclusion` are required) |
| `DB/CH_Genes.txt` | List of clonal hematopoiesis genes (column `Gene`); read by `10_summarize_fgbio.R` |
| `DB/` | CHIP gene lists and panel target regions (BED files, see below) |
| `Panel_of_normal_fgbio.xlsx` | Panel of Normals, created by `05_PanelOfNormals.R` (sheets `Background`, `Technical background` and `Noise`). **Not included:** the PoN depends on the cohort and processing conditions. Build your own from healthy donor samples processed with the same workflow |
| `PoN_Samples.xlsx` | Sample sheet for the PoN (sheet `PoN Samples`, columns `ID` and `Batch`); required only for `05_PanelOfNormals.R`. **Not included: contains sample identifiers.** |

### CHIP gene lists (`DB/`)

The CHIP classification uses three gene lists (SARAH CH gene lists, v1.2):

- `M-CHIP_Genes_modified_v1.2.tsv`
- `L-CHIP_Genes_pathogenic_modified_v1.2.tsv`
- `L-CHIP_Genes_putative_modified_v1.2.tsv`

They are stored in the `DB/` folder. Make sure the constants `DBPATH_MCHIP`, `DBPATH_LCHIP_PATH` and `DBPATH_LCHIP_PUT` at the top of `10_summarize_fgbio.R` point to these files.

Position criteria in the `CHIP_Variants` columns are interpreted as **amino acid (protein) positions** by default and are matched against `HGVSp_Short`. Only tokens with an explicit `c.` prefix are treated as cDNA positions. Position-specific criteria are only evaluated if the variant is annotated on the MANE Select transcript listed in the CHIP database (`MANE_SELECT` column of the MAF, see upstream workflow).

### Panel target regions (`DB/`)

The LION panel target regions (109 genes, hybrid-capture design) are provided as BED files for both human genome builds:

- `DB/LION_Panel_GRCh37.bed`
- `DB/LION_Panel_GRCh38.bed`

Use `DB/LION_Panel_GRCh38.bed` with the upstream workflow (`fgbio_v1_GRCh38`, variable `BEDFILE`), which is designed for GRCh38, or `DB/LION_Panel_GRCh37.bed` for data aligned to hg19/GRCh37. The panel was originally designed on hg19/GRCh37 and lifted over to GRCh38; the GRCh37 file is not used by any script in this repository. A complete gene list including annotated hotspot regions is provided in Supplementary Table S1 of the manuscript and in `DB/Gene_list.txt`.

### DrugBank: must be obtained separately

> **Important:** The DrugBank data are **not** part of this repository. DrugBank is distributed under a license that does not allow redistribution. To use the drug annotation, you must download the data yourself from <https://go.drugbank.com> (a license/account is required, depending on your use case) and place them in a local folder of your choice.

The scripts were developed with DrugBank version 5.1.12 and expect these two files:

- `all drug links.csv`
- `drugbank_all_target_polypeptide_ids/pharmacologically_active.csv`

Adjust the corresponding paths in `10_summarize_fgbio.R` (section *Data DrugBank*) to your local copy. Without these files the script will stop at that step.

---

## Usage

### 1. Run the upstream workflow

Run `fgbio_v1_GRCh38` once per sample (see [Variables](#variables)) and arrange the result folders as `../<BATCH>/fgbio_results/<SAMPLEID>/`.

### 2. Build the Panel of Normals

Requires VarDict calls of healthy donor samples (processed with `fgbio_v1_GRCh38`) and a sample sheet (`PoN_Samples.xlsx`, columns `Batch` and `ID`).

```bash
Rscript 05_PanelOfNormals.R
```

Output: `Panel_of_normal_fgbio.xlsx`

Only calls with more than two alternate reads are used. Variants are assigned to the PoN *Background* if they occur in at least two healthy donors with VAF >= 35% (constitutional background); variants occurring in at least two donors with VAF < 35% form the *Technical background*. Per-variant VAF noise thresholds (*Noise*) are calculated as `1.25 x` the highest VAF observed in the healthy donors after removing the single highest value (capped at 1; requires the variant in at least two donors).

### 3. Summarize a batch

Interactive: set `BATCHID` in `10_summarize_fgbio.R` and source the script.

Command line (overrides `BATCHID`):

```bash
Rscript 10_summarize_fgbio.R batch_a001
```

Several batches (e.g., batches 1 to 14):

```bash
for i in $(seq -w 1 14); do
    Rscript 10_summarize_fgbio.R "batch_a0${i}"
done
```

All paths in the scripts are **relative** to the working directory (`../<BATCH>/...`, `../DB/...`). Start the scripts from the directory in which they are located, or adjust the paths in the setup section of `10_summarize_fgbio.R`.

---

## Output

Written to `../<BATCH>/fgbio_tables/` and additionally copied to `../liq_bio_runs/<BATCH>/` (archive folder, set via `PATHARCHIVE`):

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

- **Filtering:** population frequency (`MAX_AF` <= 0.001 or missing), minimum alt read count (> 2), VAF above the respective reporting threshold, variant classes of interest (frameshift, in-frame indel, missense, nonsense, nonstop, splice site/region), exclusion of PoN background variants; the per-variant PoN noise threshold is annotated as `t_vaf_above_noise`. For Mutect2, only `PASS` calls are reported in the filtered table.
- **CHIP classification:** gene- and variant-specific matching against the three CHIP lists for all variants with VAF >= 1% and more than two alternate reads; the binomial test (H0: VAF = 0.5 and H0: VAF = 1.0) flags variants consistent with germline heterozygosity or homozygosity (`CHIP_binom_flag`). Known high-VAF CHIP exceptions are always called somatic.
- **Tumor fraction:** surrogate TF following Husain et al. (2022, *JCO Precis Oncol*): the maximum VAF among somatic variants (VarDict, 0.5% filtered table, above the PoN noise threshold) after exclusion of CHIP variants and putative germline variants (VAF >= 35%). Additional VAF distribution statistics (p80, mean, SD, median), an outlier flag (`TF_max_outlier`, TF_max/p80 > 3), Mutect2 PASS concordance of the max-VAF variant and an experimental KDE-based estimate (`TF_kde`) are reported. Interpret TF estimates with caution in samples with few somatic variants.
- **VAF thresholds:** the thresholds used for the filtered tables are reporting thresholds of the bioinformatic workflow. They are not validated limits of detection (see the manuscript for the analytical validation).

---

## Notes and limitations

- The PoN is built from VarDict calls of healthy donors and applied to both callers.
- The scripts were developed for the LION panel and the fgbio-based workflow; use with other panels requires adjustments of the gene list, PoN and CHIP lists.
- The code is provided for research use and has not been validated as an in-vitro diagnostic.

---

## Citation

If you use this code, please cite:

> Feierabend S, Künstner A, et al. *A liquid biopsy-centered, pan-cancer, open next generation sequencing panel to support clinical decision-making (LION panel).* [medRxiv, 2026, DOI: https://doi.org/10.64898/2026.06.05.26354976]

Code archive: this repository

The CHIP classification is adapted from the CHIP pipeline by L. Schawe (@LSchawe); the ACMG classification function is taken from the MIRACUM pipeline (https://doi.org/10.3390/cancers15133456). Please acknowledge these sources accordingly.

---

## License

This code is released under the [MIT License](LICENSE).

Third-party resources (DrugBank, CHIP gene lists) are subject to their own licenses and terms of use.

---

## Contact

Axel Künstner, University of Luebeck
