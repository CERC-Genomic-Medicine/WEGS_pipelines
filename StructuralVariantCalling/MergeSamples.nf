#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2022
*/

inputFiles = Channel.fromPath(params.inputFiles).map { file -> file.getSimpleName() }

process InitializeFileList{
    label "InitializeFileList"
    cache "lenient"
    executor "local"
    cpus 1

    input:   
    val input_label from inputFiles

    output:
    val input_label into samples

    """
    mkdir -p ${params.result}/mergeSamples
    cat /dev/null >  ${params.result}/mergeSamples/ALL_merged_FileList.list
    cat /dev/null > ${params.result}/mergeSamples/DEL_merged_FileList.list
    cat /dev/null > ${params.result}/mergeSamples/DUP_merged_FileList.list
    cat /dev/null > ${params.result}/mergeSamples/INS_merged_FileList.list
    cat /dev/null > ${params.result}/mergeSamples/INV_merged_FileList.list
    cat /dev/null > ${params.result}/mergeSamples/TRABND_merged_FileList.list
    """



}

process CreateFileList{
    label "CreateFileList"
    cache "lenient"
    executor "local"
    cpus 1

    input:   
    val sample from samples

    output:
    val sample into del, ins, dup, inv, trabnd, all

    """
    echo ${params.result}/mergeTools/${sample}/ALL_merged.vcf >> ${params.result}/mergeSamples/ALL_merged_FileList.list
    echo ${params.result}/mergeTools/${sample}/DEL_merged.vcf >> ${params.result}/mergeSamples/DEL_merged_FileList.list
    echo ${params.result}/mergeTools/${sample}/DUP_merged.vcf >> ${params.result}/mergeSamples/DUP_merged_FileList.list
    echo ${params.result}/mergeTools/${sample}/INS_merged.vcf >> ${params.result}/mergeSamples/INS_merged_FileList.list
    echo ${params.result}/mergeTools/${sample}/INV_merged.vcf >> ${params.result}/mergeSamples/INV_merged_FileList.list
    echo ${params.result}/mergeTools/${sample}/TRABND_merged.vcf >> ${params.result}/mergeSamples/TRABND_merged_FileList.list
    """
}

process Del {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "8h"

   input:   
   val(input_label) from del.collect()

   output:
   file "samples_merged_DEL.tmp.vcf" into del_merged
   
   """   
   ${params.survivor} merge ${params.result}/mergeSamples/DEL_merged_FileList.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} 0 samples_merged_DEL.tmp.vcf
   """
}

process Dup {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "8h"

   input:   
   val(input_label) from dup.collect()

   output:
   file "samples_merged_DUP.tmp.vcf" into dup_merged
   
   """   
   ${params.survivor} merge ${params.result}/mergeSamples/DUP_merged_FileList.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} 0 samples_merged_DUP.tmp.vcf
   """
}

process Ins {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "8h"

   input:   
   val(input_label) from ins.collect()

   output:
   file "samples_merged_INS.tmp.vcf" into ins_merged
   
   """   
   ${params.survivor} merge ${params.result}/mergeSamples/INS_merged_FileList.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} 0 samples_merged_INS.tmp.vcf
   """
}
process Inv {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "8h"

   input:   
   val(input_label) from inv.collect()

   output:
   file "samples_merged_INV.tmp.vcf" into inv_merged
   
   """   
   ${params.survivor} merge ${params.result}/mergeSamples/INV_merged_FileList.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} 0 samples_merged_INV.tmp.vcf
   """
}
process TraBnd {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "8h"

   input:   
   val(input_label) from trabnd.collect()

   output:
   file "samples_merged_TRABND.tmp.vcf" into trabnd_merged
   
   """   
   ${params.survivor} merge ${params.result}/mergeSamples/TRABND_merged_FileList.list ${params.breakpoint_dist} 1 ${params.use_type} ${params.use_strand} ${params.dist_based} 0 samples_merged_TRABND.tmp.vcf
   """
}

process All {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"

   input:   
   file(del) from del_merged
   file(dup) from dup_merged
   file(ins) from ins_merged
   file(inv) from inv_merged
   file(trabnd) from trabnd_merged

   output:
   tuple file("samples_merged_ALL.sorted.vcf"),file("samples_merged_DEL.vcf"),file("samples_merged_DUP.vcf"),file("samples_merged_INS.vcf"),file("samples_merged_INV.vcf"),file("samples_merged_TRABND.vcf") into merged
   
   publishDir "${params.result}/mergeSamples", mode: "copy"
   """  
   module load bcftools
   module load samtools
   module load vcftools
   module load r
   export R_LIBS=/home/praveen/R/x86_64-pc-linux-gnu-library/4.1

   cat ${del} > samples_merged_ALL.tmp.vcf
   grep -v "^#" ${dup} >> samples_merged_ALL.tmp.vcf
   grep -v "^#" ${ins} >> samples_merged_ALL.tmp.vcf
   grep -v "^#" ${inv} >> samples_merged_ALL.tmp.vcf
   grep -v "^#" ${trabnd} >> samples_merged_ALL.tmp.vcf
   
   bcftools sort -Ov samples_merged_ALL.tmp.vcf > samples_merged_ALL.sorted.vcf
   
   Rscript ${params.scripts}/renameSVIds.R samples_merged_ALL.sorted.vcf samples_merged_ALL.sorted.vcf.gz
   gunzip -f samples_merged_ALL.sorted.vcf.gz

   grep -E '^#|SVTYPE=DEL' samples_merged_ALL.sorted.vcf > samples_merged_DEL.vcf
   grep -E '^#|SVTYPE=DUP' samples_merged_ALL.sorted.vcf > samples_merged_DUP.vcf
   grep -E '^#|SVTYPE=DEL|SVTYPE=DUP' samples_merged_ALL.sorted.vcf > samples_merged_CNV.vcf
   grep -E '^#|SVTYPE=INS' samples_merged_ALL.sorted.vcf > samples_merged_INS.vcf
   grep -E '^#|SVTYPE=INV' samples_merged_ALL.sorted.vcf > samples_merged_INV.vcf
   grep -E '^#|SVTYPE=BND' samples_merged_ALL.sorted.vcf > samples_merged_BND.vcf
   sed 's/SVTYPE=BND/SVTYPE=TRA/' samples_merged_ALL.sorted.vcf | grep -E '^#|SVTYPE=TRA' > samples_merged_TRABND.vcf
   """ 
}
