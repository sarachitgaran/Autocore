<p align="center">
<img width="420" height="179" alt="image" src="https://github.com/user-attachments/assets/5ca12b17-9fbd-4e28-a141-a77d8eb6fc7d" />

# Autocore

Autocore is a bioinformatics pipeline for processing and analyzing RNA-seq data, with a focus on autoimmune diseases. The repository contains scripts and resources for data prefetching, alignment, quantification, and downstream analysis.

## Disease Modules

- **IBD (Inflammatory Bowel Disease):** Chronic inflammation of the digestive tract, including Crohn's Disease(CD) and Ulcerative Colitis(UC).
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
   - Execute the shell scripts in order (`0_prefetch.sh` → `4_htseq_count.sh`).
   - Use R scripts in disease folders for downstream analysis.

## Contributing

Contributions are welcome! Please open issues or submit pull requests for improvements or bug fixes.

## License

This project is licensed under the MIT License.

## Contact

For questions or collaboration, contact the repository owner via GitHub.
