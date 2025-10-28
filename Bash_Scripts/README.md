# Bash Scripts for Differential Gene Expression Analysis

This directory contains SLURM batch scripts used in the RNA-seq data processing pipeline for the Autocore project. SLURM is a workload manager used to schedule jobs on high-performance computing (HPC) clusters. Each script here is designed to be submitted to a SLURM-managed cluster, specifying resources such as CPUs, memory, and runtime for efficient parallel processing.

Each script corresponds to a major step in the workflow, from retreiving the data to quantification.

## Scripts Overview

- **0_prefetch.sh**: Downloads raw sequencing data using accession numbers listed in `Acc_list.txt`.
- **1_fastq_dump.sh**: Converts downloaded `.sra` files to FASTQ files using `fastq-dump`. Since the library is Paired-end, we ended up with 2 FASTQ files for each sample.
- **2_hisat.sh**: Aligns FASTQ files to the reference genome using HISAT2 aligner.
- **3_samtools.sh**: Converts SAM files to sorted and indexed BAM files using SAMtools.
- **4_htseq_count.sh**: Counts gene features from sorted BAM files using HTSeq.

## Usage

Run each script in order as your data progresses through the pipeline.

To submit a script to SLURM, use:
```sh
sbatch script_name.sh
```

---

For questions or issues, please refer to the main project README or contact the repository owner.
