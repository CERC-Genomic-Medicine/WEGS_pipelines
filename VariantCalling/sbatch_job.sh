#!/bin/bash
#SBATCH --job-name=vc_gatk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00

module load nextflow
module load singularity
module load bcftools
nextflow run /home/praveen/projects/def-vmooser/praveen/ExomePlus/VariantCalling/Pipeline.nf -w /home/praveen/scratch/wes+wgs -resume 6354c42a-335a-4bcd-baba-141c5f750f18