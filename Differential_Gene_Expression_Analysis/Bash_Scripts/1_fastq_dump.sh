#!/bin/bash
#SBATCH --job-name=fastq_dump
#SBATCH --cpus-per-task=1
#SBATCH --time=01-00:00:00
#SBATCH --mem=250GB
#SBATCH --output=fastq_dump_%j.out
#SBATCH --error=fastq_dump_%j.err
#SBATCH --mail-type=END,FAIL,TIME_lIMIT
#SBATCH --mail-user=sara@lbi-netmed.com

module load sratoolkit
for i in $(seq 33775194 33775216); do
	id="SRR${i}"
    dir="$id"
    sra_file="$dir/$id.sra"
    echo "Processing $sra_file..."
    fastq-dump --split-files --gzip --outdir "$dir" "$sra_file"
done

# Alternative reproducible approach: loop over all SRR*/SRR*.sra files automatically
# This makes the script independent of specific numeric ranges and more robust.
for sra_file in SRR*/SRR*.sra; do
    id=$(basename "$sra_file" .sra)
    dir=$(dirname "$sra_file")
    echo "Processing $sra_file..."
    fastq-dump --split-files --gzip --outdir "$dir" "$sra_file"
done

