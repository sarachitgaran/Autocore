# Autocore

Autocore is a bioinformatics pipeline for processing and analyzing RNA-seq data, with a focus on autoimmune diseases. The repository contains scripts and resources for data prefetching, alignment, quantification, and downstream analysis.

## Repository Structure

- `0_prefetch.sh` — Script for prefetching raw sequencing data.
- `1_fastq_dump.sh` — Converts SRA files to FASTQ format.
- `2_hisat.sh` — Aligns reads using HISAT2.
- `3_samtools.sh` — Processes alignment files with SAMtools.
- `4_htseq_count.sh` — Quantifies gene expression using HTSeq.
- `counts/` — Contains gene count files for various samples.
- `IBD/`, `Psoriasis/`, `RA/`, `Sjörgen/`, `SLE/`, `T1D/` — Disease-specific folders with metadata, counts, plots, and analysis scripts.
- `README.md` — Project documentation.
- `Autocore.Rproj` — R project file for reproducible analysis.
- `renv/` — R environment management.

## Getting Started

1. **Clone the repository:**
   ```sh
   git clone https://github.com/sarachitgaran/Autocore.git
   cd Autocore
   ```
2. **Install dependencies:**
   - R (recommended: use RStudio)
   - HISAT2, SAMtools, HTSeq
   - [renv](https://rstudio.github.io/renv/) for R package management

3. **Run the pipeline:**
   - Execute the shell scripts in order (`0_prefetch.sh` → `4_htseq_count.sh`).
   - Use R scripts in disease folders for downstream analysis.

## Data Organization

- **counts/**: Contains raw and processed count files for each sample.
- **Disease folders**: Each folder contains metadata, counts, enrichment tests, plots, and tables for a specific autoimmune disease.

## Contributing

Contributions are welcome! Please open issues or submit pull requests for improvements or bug fixes.

## License

This project is licensed under the MIT License.

## Contact

For questions or collaboration, contact the repository owner via GitHub.
