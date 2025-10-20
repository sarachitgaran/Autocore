#!/bin/bash
#SBATCH --job-name=samtools_convert
#SBATCH --cpus-per-task=40
#SBATCH --time=05-00:00:00
#SBATCH --mem=500GB
#SBATCH --output=samtools_conversion_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_LIMIT_50,TIME_lIMIT
#SBATCH --mail-user=sara@lbi-netmed.com

module load samtools

for id in $(seq 1721280 1721313); do
    sample="SRR${id}"
    echo "Processing $sample ..."
    samtools view -@ 40 -bS "./${sample}.sam" | \
    samtools sort -@ 40 -o "./${sample}_sorted.bam"
    samtools index "./${sample}_sorted.bam"
done
