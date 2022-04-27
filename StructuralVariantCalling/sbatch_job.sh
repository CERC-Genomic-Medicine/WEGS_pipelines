#!/bin/bash
#SBATCH --job-name=Lumpy
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=4096
#SBATCH --cpus-per-task=1
#SBATCH --time=50:00:00

module load nextflow
nextflow run /home/praveen/projects/rrg-vmooser/praveen/pipelines/ExomePlus/StructuralVariantCalling/Pipeline.nf -resume 664b2fb8-6e4c-4d47-9fb6-cbf84662d71a