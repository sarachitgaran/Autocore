#!/bin/bash
#
#SBATCH --job-name=samtools_convert
#SBATCH --cpus-per-task=40
#SBATCH --time=03-00:00:00
#SBATCH --mem=250GB
#SBATCH --output=samtools_conversion_%j.out
#SBATCH --mail-type=END,FAIL,TIME_LIMIT_50,TIME_lIMIT
#SBATCH --mail-user=sara@lbi-netmed.com

module load HISAT2
module load samtools
mkdir hisat2_output 

for id in $(seq 33775194 33775216); do
    sample="SRR${id}"
    hisat2 -x /lisc/scratch/menche/sara/grch38/genome \
        -1 "/lisc/scratch/menche/sara/SLE/${sample}/${sample}_1.fastq.gz" \
        -2 "/lisc/scratch/menche/sara/SLE/${sample}/${sample}_2.fastq.gz" \
        -S "hisat2_output/${sample}.sam"
done



# Alternative reproducible approach: loop over all paired-end FASTQ files in SRR* directories
# This makes the script independent of specific numeric ranges and more robust.
for dir in SRR*/; do
    sample=$(basename "$dir")
    fq1="$dir/${sample}_1.fastq.gz"
    fq2="$dir/${sample}_2.fastq.gz"
    if [[ -f "$fq1" && -f "$fq2" ]]; then
        echo "Processing $sample ..."
        hisat2 -x /lisc/scratch/menche/sara/grch38/genome \
            -1 "$fq1" \
            -2 "$fq2" \
            -S "hisat2_output/${sample}.sam"
    fi
done


