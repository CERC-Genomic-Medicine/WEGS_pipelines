#!/bin/bash
#SBATCH --job-name=vc_gatk
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00:00

module load nextflow
module load singularity
module load bcftools
module load samtools
nextflow run /home/praveen/projects/def-vmooser/praveen/ExomePlus/VariantCalling/Pipeline.nf -w /home/praveen/scratch/wes+wgs -resume eac1361a-05fe-467a-b6f0-87d0fe78dd37