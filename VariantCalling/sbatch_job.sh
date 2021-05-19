#!/bin/bash
#SBATCH --job-name=vc_gatk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00

module load nextflow
module load singularity
module load bcftools
module load samtools
nextflow run /home/praveen/projects/def-vmooser/praveen/ExomePlus/VariantCalling/NoRecalPipeline.nf -w /home/praveen/scratch/wes+wgs