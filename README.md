# Autocore

Autocore is a comprehensive bioinformatics pipeline designed for the processing and analysis of RNA-seq data, with an emphasis on autoimmune diseases. The repository provides modular scripts and resources for data prefetching, alignment, quantification, and downstream analysis.

---

## Disease Modules

Autocore supports the following autoimmune diseases:

- **IBD (Inflammatory Bowel Disease):** Chronic inflammation of the digestive tract, including Crohn's disease and ulcerative colitis.
- **Psoriasis:** Autoimmune skin disorder causing red, scaly patches.
- **RA (Rheumatoid Arthritis):** Autoimmune disorder affecting joints, causing pain and swelling.
- **Sjörgen's Syndrome:** Autoimmune disease targeting moisture-producing glands, leading to dry mouth and eyes.
- **SLE (Systemic Lupus Erythematosus):** Multi-system autoimmune disease affecting skin, joints, kidneys, and other organs.

Each disease module includes:
- Disease-specific metadata
- Count files
- Analysis scripts
- Plots and tables

---

## Getting Started

### 1. Clone the Repository

```sh
git clone https://github.com/sarachitgaran/Autocore.git
cd Autocore
```

### 2. Install Dependencies

- **R** (recommended: use <a href="https://posit.co/download/rstudio-desktop/">RStudio</a>)
- **HISAT2**, **SAMtools**, **HTSeq**
- <a href="https://rstudio.github.io/renv/">renv</a> for R package management

### 3. Run the Pipeline

- Execute the shell scripts sequentially, e.g.:
  - `0_prefetch.sh` → `4_htseq_count.sh`
- Use the R scripts within each disease module folder for downstream analysis.

---

## Contributing

We welcome contributions! If you would like to report issues or propose improvements, please open an issue or submit a pull request.

---

## License

This project is licensed under the MIT License.

---

## Contact

For questions, feedback, or collaboration inquiries, please contact the repository owner via <a href="https://github.com/sarachitgaran/Autocore/issues">GitHub Issues</a> or by opening a discussion.
