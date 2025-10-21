# Autocore

Autocore is a bioinformatics pipeline for processing and analyzing RNA-seq data, with a focus on autoimmune diseases. The repository contains scripts and resources for data prefetching, alignment, quantification, and downstream analysis.

## Repository Structure (Recommended)

```
Autocore/
│
├── data/                # Raw and processed data (input/output)
│   ├── raw/
│   └── processed/
│
├── scripts/             # All pipeline scripts (shell, R, Python, etc.)
│   ├── prefetch/
│   ├── fastq/
│   ├── alignment/
│   ├── quantification/
│   └── analysis/
│
├── results/             # Output results (counts, plots, tables)
│   ├── counts/
│   ├── plots/
│   └── tables/
│
├── metadata/            # Metadata files (e.g., sample info, experimental design)
│
├── docs/                # Documentation, manuals, and references
│
├── env/                 # Environment and config files (renv, conda, etc.)
│
├── disease_modules/     # Disease-specific analyses
│   ├── IBD/
│   ├── Psoriasis/
│   ├── RA/
│   ├── Sjörgen/
│   └── SLE/
│
├── README.md
├── LICENSE
└── Autocore.Rproj
```

## Disease Modules

- **IBD (Inflammatory Bowel Disease):** Chronic inflammation of the digestive tract, including Crohn's disease and ulcerative colitis.
- **Psoriasis:** Autoimmune skin disorder causing red, scaly patches.
- **RA (Rheumatoid Arthritis):** Autoimmune disorder affecting joints, causing pain and swelling.
- **Sjörgen's Syndrome:** Autoimmune disease targeting moisture-producing glands, leading to dry mouth and eyes.
- **SLE (Systemic Lupus Erythematosus):** Multi-system autoimmune disease affecting skin, joints, kidneys, and other organs.

Each disease module contains:
- Disease-specific metadata
- Count files
- Analysis scripts
- Plots and tables

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
   - Organize your files according to the structure above.
   - Execute scripts in the `scripts/` folder stepwise.
   - Use disease modules for downstream analysis and visualization.

## Contributing

Contributions are welcome! Please open issues or submit pull requests for improvements or bug fixes.

## License

This project is licensed under the MIT License.

## Contact

For questions or collaboration, contact the repository owner via GitHub.
