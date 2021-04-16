#!/bin/bash
#SBATCH --job-name=merge_bams
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=16:00:00

module load nextflow
module load samtools
module load picard
module load bedtools
nextflow run /home/praveen/projects/def-vmooser/praveen/ExomePlus/WES+WGS/MergeBAMsPipeline.nf -w /home/praveen/scratch/wes+wgs