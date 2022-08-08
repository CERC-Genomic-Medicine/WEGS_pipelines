#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2021
*/

process BaseRecalibrator {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "24h"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref -B ${params.bundle}:/bundle"

   input:   
   tuple val(input_label), file(input), file(index) from Channel.fromPath(params.inputFilesToRecal).map{ bam -> [ bam.getSimpleName(), bam, bam + (bam.getExtension() == "bam" ? ".bai" : ".crai") ] }

   output:
   tuple val(input_label), file(input), file(index), file("${input_label}.table") into for_ApplyBQSR



   """
   gatk --java-options "-Xmx14G" BaseRecalibrator -I ${input} -R /ref/${params.referenceGenome}.fa --known-sites /bundle/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /bundle/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O ${input_label}.table
   """
}


process ApplyBQSR {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "24h"
   scratch '$SLURM_TMPDIR'

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref -B ${params.recalFolder}:/recal"

   input:   
   tuple val(input_label), file(input), file(index), file(recalTable) from for_ApplyBQSR
   
   script:
   if ( params.inputFileType == "bam" )
      """
      gatk --java-options "-Xmx14G" ApplyBQSR -I ${input} -R /ref/${params.referenceGenome}.fa --bqsr-recal-file ${recalTable} -O ${input_label}.recalV4.bam
      """
   else if ( params.inputFileType == "cram" )
      """
      gatk --java-options "-Xmx14G" ApplyBQSR -I ${input} -R /ref/${params.referenceGenome}.fa --bqsr-recal-file ${recalTable} -O ${input_label}.recalV4.cram
      rm ${input_label}.recalV4.*bai
      samtools index ${input_label}.recalV4.cram
      cp ${input_label}.recalV4.* /recal
      """
   else
      throw new IllegalArgumentException("Unknown input file type $params.inputFileType") 
}