#!/bin/bash
#SBATCH --job-name=prefetch
#SBATCH --cpus-per-task=1
#SBATCH --time=05-00:00:00
#SBATCH --mem=250GB
#SBATCH --output=prefetch_%j.out
#SBATCH --error=prefetch_%j.err
#SBATCH --mail-type=END,FAIL,TIME_lIMIT
#SBATCH --mail-user=sara@lbi-netmed.com

module load sratoolkit

prefetch --option-file Acc_list.txt
