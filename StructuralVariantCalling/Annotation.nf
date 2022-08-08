#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2022
*/

process Annotate {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "6h"

   """
   module load r
   module load samtools
   module load bcftools
   module load bedtools

   export R_LIBS=/home/praveen/R/x86_64-pc-linux-gnu-library/4.1
   export ANNOTSV=/home/praveen/projects/rrg-vmooser/praveen/tools/SV/AnnotSV/AnnotSV
   
   mkdir -p ${params.result}/AnnotSV
   bgzip -d -c ${params.result}/filtering/fromSMOOVE/samples_merged_ALL.Final.vcf.gz > samples_merged_ALL.Final.vcf

   ${params.survivor} vcftobed samples_merged_ALL.Final.vcf -99999999 99999999 samples_merged_ALL.Final.bed
   cat samples_merged_ALL.Final.bed | cut -f1,2,5,7,11 | awk -F '\t' '{ \$3 = (\$3 == "0" ? \$2+1 : \$3) } 1' OFS='\t'| awk '(\$3 > \$2 )' > samples_merged_ALL.Final.clean.bed

   ${params.annotsv} -SVinputFile samples_merged_ALL.Final.clean.bed -outputDir ${params.result}/AnnotSV -SVinputInfo 1 -reciprocal 1 -svtBEDcol 5

   Rscript ${params.scripts}/complementAnnotation.R ${params.result}/AnnotSV/samples_merged_ALL.Final.clean.annotated.tsv ${params.result}/AnnotSV/samples_merged_ALL.Final.clean.annotated.plus.tsv h38
   """
}