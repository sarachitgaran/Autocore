#!/bin/bash
#
#SBATCH --job-name=htseq_count
#SBATCH --cpus-per-task=1
#SBATCH --time=02-00:00:00
#SBATCH --mem=128GB
#SBATCH --output=htseq_count_%j.out
#SBATCH --mail-type=BEGIN,END,FAIL,TIME_lIMIT
#SBATCH --mail-user=sara@lbi-netmed.com

module load htseq

for bamfile in *_sorted.bam; do
    sample=${bamfile%%_sorted.bam}
    htseq-count \
        --format=bam \
        --order=pos \
        --stranded=no \
        --type=exon \
        --idattr=gene_id \
        "$bamfile" \
        *.gtf \
        > "${sample}_counts.txt"
done

# Alternative reproducible approach: loop over all *_sorted.bam files for counting
# This makes the script independent of hardcoded values and more robust.
for bamfile in SRR*_sorted.bam; do
    sample=${bamfile%%_sorted.bam}
    htseq-count \
        --format=bam \
        --order=pos \
        --stranded=no \
        --type=exon \
        --idattr=gene_id \
        "$bamfile" \
        *.gtf \
        > "${sample}_counts.txt"
done
