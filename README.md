<p align="center">
<img width="420" height="179" alt="image" src="https://github.com/user-attachments/assets/5ca12b17-9fbd-4e28-a141-a77d8eb6fc7d" />

# Autocore

Autocore is a bioinformatics pipeline for processing and analyzing RNA-seq data, with a focus on autoimmune diseases. The repository contains scripts and resources for data prefetching, alignment, quantification, and downstream analysis.

## Disease Modules

- **IBD (Inflammatory Bowel Disease):** Chronic inflammation of the digestive tract, including Crohn's Disease(CD) and Ulcerative Colitis(UC) which can ultimatley lead to Colorectal Cancer (CRC).
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
   - R
   - HISAT2, SAMtools, HTSeq

3. **Run the pipeline:**
   - Execute the shell scripts in order (`0_prefetch.sh` → `4_htseq_count.sh`).
   - Use R scripts in disease folders for downstream analysis.

## Repository Structure

All Linux-based processing scripts are located in the `Bash_Scripts` directory, with each step of the pipeline implemented as a separate SLURM batch file. The resulting gene count files are collected in the  `Count_Files` directory. Differential gene expression (DGE) analysis and its results—including plots and tables—are organized in the `Differential_Gene_Expression_Analysis` directory. Further downstream analyses are stored in the `Downstream_Analysis` drectory.

## Contributing

Contributions are welcome! Please open issues or submit pull requests for improvements or bug fixes.
