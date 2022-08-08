#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2022
*/

inputFiles = Channel.fromPath(params.inputFiles).map { file -> file.getSimpleName() }

process MakeWorkDir{
    errorStrategy 'retry'
    maxRetries 3
    cache "lenient"
    cpus 1
    memory "4GB"
    time "1h"

    input:   
    val input_label from inputFiles

    output:
    val input_label into del, ins, dup, inv, trabnd

    """
    mkdir -p ${params.result}/mergeTools/${input_label}
    """



}

process Del {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "4h"

   input:   
   val(input_label) from del

   output:
   val input_label into del_merged
   
   """
   module load bcftools
   module load samtools
   module load vcftools
   module load r
   
   export R_LIBS=/home/praveen/R/x86_64-pc-linux-gnu-library/4.1
   find ${params.result}/prepareVCFsToMerge/ -type f -name "${input_label}.manta.vcf.filt.DEL" > DEL_filtered.vcf.list1
   ${params.survivor} merge DEL_filtered.vcf.list1 ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} DEL_merged.vcf1
   sed -i 's|FORMAT=<ID=DR,Number=1,Type=Integer|FORMAT=<ID=DR,Number=1,Type=String|g' DEL_merged.vcf1
   sed -i 's|ID=LN,Number=1,Type=Integer|ID=LN,Number=1,Type=String|g' DEL_merged.vcf1

   find ${params.result}/prepareVCFsToMerge/ -type f -name "${input_label}.*.vcf.filt.DEL" > DEL_filtered.vcf.list2
   ${params.survivor} merge DEL_filtered.vcf.list2 ${params.breakpoint_dist} 2 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} DEL_merged.vcf2

   find . -type f -name "DEL_merged.vcf[1|2]" > DEL_filtered.vcf.list
   ${params.survivor} merge DEL_filtered.vcf.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} DEL_merged.vcf3

   vcf-sort -c DEL_merged.vcf3 > DEL_merged.sorted.vcf3
   bgzip -f -c DEL_merged.sorted.vcf3 > DEL_merged.sorted.vcf3.gz
   tabix -p vcf DEL_merged.sorted.vcf3.gz
   Rscript ${params.scripts}/fixSURVIVORgenotypes.R DEL_merged.sorted.vcf3.gz ${input_label} DEL_merged.vcf.gz
   gunzip -f DEL_merged.vcf.gz
   cp DEL_merged.vcf ${params.result}/mergeTools/${input_label}
   """
}

process Ins {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "4h"

   input:   
   val(input_label) from ins

   output:
   val input_label into ins_merged
   
   """
   module load bcftools
   module load samtools
   module load vcftools
   module load r
   
   export R_LIBS=/home/praveen/R/x86_64-pc-linux-gnu-library/4.1
   find ${params.result}/prepareVCFsToMerge/ -type f -name "${input_label}.manta.vcf.filt.INS"  > INS_filtered.vcf.list
   find ${params.result}/prepareVCFsToMerge/ -type f -name "${input_label}.breakseq.vcf.filt.INS"  >> INS_filtered.vcf.list
   ${params.survivor} merge INS_filtered.vcf.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} INS_merged.vcf1
   sed -i 's|FORMAT=<ID=DR,Number=1,Type=Integer|FORMAT=<ID=DR,Number=1,Type=String|g' INS_merged.vcf1
   sed -i 's|ID=LN,Number=1,Type=Integer|ID=LN,Number=1,Type=String|g' INS_merged.vcf1

   vcf-sort -c INS_merged.vcf1 > INS_merged.sorted.vcf1
   bgzip -f -c INS_merged.sorted.vcf1 > INS_merged.sorted.vcf1.gz
   tabix -p vcf INS_merged.sorted.vcf1.gz
   Rscript ${params.scripts}/fixSURVIVORgenotypes.R INS_merged.sorted.vcf1.gz ${input_label} INS_merged.vcf.gz
   gunzip -f INS_merged.vcf.gz
   cp INS_merged.vcf ${params.result}/mergeTools/${input_label}
   """
}

process Dup {
   errorStrategy 'retry'
   maxRetries 3
   maxRetries 1
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "4h"

   input:   
   val(input_label) from dup

   output:
   val input_label into dup_merged
   
   """
   module load bcftools
   module load samtools
   module load vcftools
   module load r
   
   export R_LIBS=/home/praveen/R/x86_64-pc-linux-gnu-library/4.1
   find ${params.result}/prepareVCFsToMerge/ -type f -name "${input_label}.*.vcf.filt.DUP" > DUP_filtered.vcf.list
   ${params.survivor} merge DUP_filtered.vcf.list ${params.breakpoint_dist} 2 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} DUP_merged.vcf1

   vcf-sort -c DUP_merged.vcf1 > DUP_merged.sorted.vcf1
   bgzip -f -c DUP_merged.sorted.vcf1 > DUP_merged.sorted.vcf1.gz
   tabix -p vcf DUP_merged.sorted.vcf1.gz
   Rscript ${params.scripts}/fixSURVIVORgenotypes.R DUP_merged.sorted.vcf1.gz ${input_label} DUP_merged.vcf.gz
   gunzip -f DUP_merged.vcf.gz
   cp DUP_merged.vcf ${params.result}/mergeTools/${input_label}
   """
}

process Inv {
   errorStrategy 'retry'
   maxRetries 3
   maxRetries 1
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "4h"

   input:   
   val(input_label) from inv

   output:
   val input_label into inv_merged
   
   """
   module load bcftools
   module load samtools
   module load vcftools
   module load r
   
   export R_LIBS=/home/praveen/R/x86_64-pc-linux-gnu-library/4.1
   find ${params.result}/prepareVCFsToMerge/ -type f -name "${input_label}.*.vcf.filt.INV" > INV_filtered.vcf.list
   ${params.survivor} merge INV_filtered.vcf.list ${params.breakpoint_dist} 2 ${params.use_type} ${params.use_strand} ${params.dist_based} ${params.min_sv_size} INV_merged.vcf1

   vcf-sort -c INV_merged.vcf1 > INV_merged.sorted.vcf1
   bgzip -f -c INV_merged.sorted.vcf1 > INV_merged.sorted.vcf1.gz
   tabix -p vcf INV_merged.sorted.vcf1.gz
   Rscript ${params.scripts}/fixSURVIVORgenotypes.R INV_merged.sorted.vcf1.gz ${input_label} INV_merged.vcf.gz
   gunzip -f INV_merged.vcf.gz
   cp INV_merged.vcf ${params.result}/mergeTools/${input_label}
   """
}

process TraBnd {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "4h"

   input:   
   val(input_label) from trabnd

   output:
   val input_label into trabnd_merged
   
   """
   module load bcftools
   module load samtools
   module load vcftools
   module load r
   module load htslib/1.11
   export R_LIBS=/home/praveen/R/x86_64-pc-linux-gnu-library/4.1
   find ${params.result}/prepareVCFsToMerge/ -type f -name "${input_label}.*.vcf.filt.BND" > TRABND_filtered.vcf.list
   find ${params.result}/prepareVCFsToMerge/ -type f -name "${input_label}.*.vcf.filt.TRA" >> TRABND_filtered.vcf.list
   ${params.survivor} merge TRABND_filtered.vcf.list 1000 1 0 1 0 0 TRABND_merged.vcf1

   vcf-sort -c TRABND_merged.vcf1 > TRABND_merged.sorted.vcf1
   bgzip -f -c TRABND_merged.sorted.vcf1 > TRABND_merged.sorted.vcf1.gz
   tabix -p vcf TRABND_merged.sorted.vcf1.gz
   Rscript ${params.scripts}/fixSURVIVORgenotypes.R TRABND_merged.sorted.vcf1.gz ${input_label} TRABND_merged.vcf.gz
   gunzip -f TRABND_merged.vcf.gz
   cp TRABND_merged.vcf ${params.result}/mergeTools/${input_label}
   """
}


process All {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "8h"

   input:   
   val(input_label) from del_merged
   val(input_label1) from ins_merged
   val(input_label2) from dup_merged
   val(input_label3) from inv_merged
   val(input_label4) from trabnd_merged

   
   """
   module load bcftools
   module load samtools
   module load vcftools
   module load htslib/1.11

   bgzip -f -c ${params.result}/mergeTools/${input_label}/DEL_merged.vcf > DEL_merged.vcf.gz
   tabix -p vcf DEL_merged.vcf.gz
   bgzip -f -c ${params.result}/mergeTools/${input_label}/INS_merged.vcf > INS_merged.vcf.gz
   tabix -p vcf INS_merged.vcf.gz
   bgzip -f -c ${params.result}/mergeTools/${input_label}/DUP_merged.vcf > DUP_merged.vcf.gz
   tabix -p vcf DUP_merged.vcf.gz
   bgzip -f -c ${params.result}/mergeTools/${input_label}/INV_merged.vcf > INV_merged.vcf.gz
   tabix -p vcf INV_merged.vcf.gz
   bgzip -f -c ${params.result}/mergeTools/${input_label}/TRABND_merged.vcf > TRABND_merged.vcf.gz
   tabix -p vcf TRABND_merged.vcf.gz  

   bcftools concat -a -O z -o ALL_merged.vcf.gz DEL_merged.vcf.gz INS_merged.vcf.gz DUP_merged.vcf.gz INV_merged.vcf.gz TRABND_merged.vcf.gz
   tabix -p vcf  ALL_merged.vcf.gz
   zcat ALL_merged.vcf.gz > ALL_merged.vcf
   cp ALL_merged.vcf ${params.result}/mergeTools/${input_label}
   """
}
