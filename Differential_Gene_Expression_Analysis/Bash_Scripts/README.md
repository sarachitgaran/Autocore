# Bash Scripts for Differential Gene Expression Analysis

This directory contains bash scripts used in the RNA-seq data processing pipeline for the Autocore project. Each script corresponds to a major step in the workflow, from data prefetching to quantification.

## Scripts Overview

- **0_prefetch.sh**: Downloads raw sequencing data using accession numbers listed in `Acc_list.txt`.
- **1_fastq_dump.sh**: Converts downloaded `.sra` files to paired-end compressed FASTQ files using `fastq-dump`. Includes both a hardcoded loop and a reproducible loop that automatically processes all available SRR directories.
- **2_hisat.sh**: Aligns paired-end FASTQ files to the reference genome using HISAT2. Contains both a hardcoded loop for specific samples and a reproducible loop for all detected SRR directories.
- **3_samtools.sh**: Converts SAM files to sorted and indexed BAM files using SAMtools. Includes both a hardcoded loop and a reproducible loop for all `.sam` files.
- **4_htseq_count.sh**: Counts gene features from sorted BAM files using HTSeq. Contains both a hardcoded loop and a reproducible loop for all sorted BAM files.

## Reproducibility

Each script now includes two approaches:
- The original approach, which uses hardcoded sample numbers for looping.
- A reproducible approach, which automatically detects and processes all relevant files or directories, making the workflow robust and adaptable to new datasets.

## Usage

Run each script in order as your data progresses through the pipeline. You may use either the original or reproducible loop, but the reproducible approach is recommended for most use cases.

---

For questions or issues, please refer to the main project README or contact the repository owner.
