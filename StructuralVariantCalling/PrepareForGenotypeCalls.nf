#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2022
*/

process PreProcessing{
    label "PreProcessing"
    cache "lenient"
    memory "8 GB"
    time "4h"
    cpus 1

    """
    module load r

    mkdir -p ${params.result}/genotype
    ${params.survivor} filter ${params.result}/mergeSamples/samples_merged_ALL.sorted.vcf NA 50 50000000 0 -1 samples_merged_ALL.filt.vcf
    grep -E '^#|SVTYPE=DEL' samples_merged_ALL.filt.vcf | cut -f1,2,3,4,5,6,7,8,9 > samples_merged_DEL.vcf
    grep -E '^#|SVTYPE=DUP' samples_merged_ALL.filt.vcf | cut -f1,2,3,4,5,6,7,8,9 > samples_merged_DUP.vcf
    grep -E '^#|SVTYPE=INV' samples_merged_ALL.filt.vcf | cut -f1,2,3,4,5,6,7,8,9 > samples_merged_INV.vcf
    
    Rscript ${params.scripts}/fixCIPOS2.R samples_merged_DEL.vcf ${params.result}/genotype/samples_merged_DEL.vcf.gz
    gunzip -f ${params.result}/genotype/samples_merged_DEL.vcf.gz
    mkdir -p ${params.result}/genotype/DELs

    Rscript ${params.scripts}/fixCIPOS2.R samples_merged_DUP.vcf ${params.result}/genotype/samples_merged_DUP.vcf.gz
    gunzip -f ${params.result}/genotype/samples_merged_DUP.vcf.gz
    mkdir -p ${params.result}/genotype/DUPs

    Rscript ${params.scripts}/fixCIPOS2.R samples_merged_INV.vcf ${params.result}/genotype/samples_merged_INV.vcf.gz
    gunzip -f ${params.result}/genotype/samples_merged_INV.vcf.gz
    mkdir -p ${params.result}/genotype/INVs
    """
}