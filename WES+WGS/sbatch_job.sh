#!/bin/bash
#SBATCH --job-name=merge_bams
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00

module load nextflow
module load samtools
module load bedtools
nextflow run /path/to/Pipeline.nf -w /path/to/working/directory