#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2021
*/

inputFiles1 = Channel.empty()
inputFiles2 = Channel.empty()

if( params.inputFileType == "bam" ) {
inputFiles1 = Channel.fromPath(params.inputFiles).map { file -> [file.getSimpleName(), file, file.getBaseName() + ".bai"] }
inputFiles2 = Channel.fromPath(params.inputFiles).map { file -> [file.getSimpleName(), file, file.getBaseName() + ".bai"] }
}
else if( params.inputFileType == "cram" ) {
inputFiles1 = Channel.fromPath(params.inputFiles).map { file -> [file.getSimpleName(), file, file + ".crai"] }
inputFiles2 = Channel.fromPath(params.inputFiles).map { file -> [file.getSimpleName(), file, file + ".crai"] }
}

intervals1 = Channel.fromPath(params.intervals).map { file -> [file.getParent().toString().split('/').last(), file]}
intervals2 = Channel.fromPath(params.intervals).map { file -> [file.getParent().toString().split('/').last(), file]}

process BaseRecalibrator {
   errorStrategy "finish"
   cache "lenient"
   cpus 1
   memory "5 GB"
   time "4h"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref -B ${params.bundle}:/bundle"

   input:   
   tuple val(input_label), file(input), file(index) from inputFiles1

   output:
   tuple val(input_label), file(input), file(index), file("${input_label}.table") into for_ApplyBQSR

   when:
   params.doBQSR

   """
   gatk --java-options "-Xmx4G" BaseRecalibrator -I ${input} -R /ref/${params.referenceGenome} --known-sites /bundle/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /bundle/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -O ${input_label}.table
   """
}

process ApplyBQSR {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "5 GB"
   time "8h"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref -B ${params.resultFolder}:/result"

   input:   
   tuple val(input_label), file(input), file(index), file(recalTable) from for_ApplyBQSR

   output:
   tuple val(input_label), file("${input_label}.recal.${params.inputFileType}") into for_RecalibratedHaplotypeCaller

   publishDir "${params.resultFolder}/recal", pattern: "${input_label}.recal.${params.inputFileType}*", mode: "copy"
   
   when:
   params.doBQSR
   
   script:
   if ( params.inputFileType == "bam" )
      """
      mkdir -p /result/recal
      gatk --java-options "-Xmx4G" ApplyBQSR -I ${input} -R /ref/${params.referenceGenome} --bqsr-recal-file ${recalTable} -O ${input_label}.recal.bam
      """
   else if ( params.inputFileType == "cram" )
      """
      mkdir -p /result/recal
      gatk --java-options "-Xmx4G" ApplyBQSR -I ${input} -R /ref/${params.referenceGenome} --bqsr-recal-file ${recalTable} -O ${input_label}.recal.cram
      """
   else
      throw new IllegalArgumentException("Unknown input file type $params.inputFileType") 
}

process RecalibratedHaplotypeCaller {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "5 GB"
   time "4h"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref"

   input: 
   tuple val(interval_label), file(interval_list), val(input_label), file(input) from intervals1.combine(for_RecalibratedHaplotypeCaller)

   output:
   tuple val(input_label), file("${interval_label}.${input_label}.vcf.gz"), file("${interval_label}.${input_label}.vcf.gz.tbi") into interval_vcfs

   when:
   params.doBQSR

   """
   gatk --java-options -Xmx4G HaplotypeCaller --native-pair-hmm-threads 1 -G StandardAnnotation -G StandardHCAnnotation -L ${interval_list} -R /ref/${params.referenceGenome} -I ${input} -O ${interval_label}.${input_label}.vcf.gz
   """

}


process HaplotypeCaller {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "8h"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref"

   input:   
   tuple val(interval_label), file(interval_list), val(input_label), file(input), file(index) from intervals2.combine(inputFiles2)

   output:
   tuple val(input_label), file("${interval_label}.${input_label}.vcf.gz"), file("${interval_label}.${input_label}.vcf.gz.tbi") into interval_vcfs2

   when:
   !params.doBQSR

   """
   gatk --java-options -Xmx12G HaplotypeCaller --native-pair-hmm-threads 1 -G StandardAnnotation -G StandardHCAnnotation -L ${interval_list} -R /ref/${params.referenceGenome} -I ${input} -O ${interval_label}.${input_label}.vcf.gz
   """
}

process RecalibratedMerge {
   errorStrategy "finish"
   cache "lenient"
   cpus 1
   memory "4 GB"
   time "3h"

   input:
   tuple val(label), file(vcfs), file(vcf_indices) from interval_vcfs.groupTuple(by: 0)

   output:
   tuple val(label), file("${label}.vcf.gz"), file("${label}.vcf.gz.tbi") into vcfs

   //publishDir "${params.resultFolder}", pattern: "${label}.vcf.gz*"
   
   when:
   params.doBQSR

   """
   find . -name "*.vcf.gz" | sort > files.txt
   bcftools concat -f files.txt -Oz -o ${label}.vcf.gz
   tabix ${label}.vcf.gz
   """
}

process Merge {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "6h"

   input:
   tuple val(label), file(vcfs), file(vcf_indices) from interval_vcfs2.groupTuple(by: 0)

   output:
   tuple val(label), file("${label}.vcf.gz"), file("${label}.vcf.gz.tbi") into vcfs2

   //publishDir "${params.resultFolder}", pattern: "${label}.vcf.gz*"

   """
   find . -name "*.vcf.gz" | sort > files.txt
   bcftools concat -f files.txt -Oz -o ${label}.vcf.gz
   tabix ${label}.vcf.gz
   """
}

process RecalibratedFilter {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "4 GB"
   time "1h"

   input:
   tuple val(label),file(vcf), file(vcf_index) from vcfs

   output:
   tuple val(label),file("${label}.hard_filter.vcf.gz"), file("${label}.hard_filter.vcf.gz.tbi") into hard_filtered_vcfs

   publishDir "${params.resultFolder}/vcfs", pattern: "${label}.hard_filter.vcf.gz*", mode: "move"

   when:
   params.doBQSR

   """
   mkdir -p ${params.resultFolder}/vcfs
   bcftools filter -s FAIL -i '(TYPE="snp" & INFO/QD>=2 & INFO/FS<=60 & INFO/MQ>=40 & INFO/SOR<=3 & (INFO/MQRankSum="." | INFO/MQRankSum>=-12.4) & (INFO/ReadPosRankSum="." | INFO/ReadPosRankSum>=-8.0)) | (TYPE="indel" & INFO/QD>=2 & INFO/FS<=200 & INFO/SOR<=10 & (INFO/ReadPosRankSum="." | INFO/ReadPosRankSum>=-20))' ${vcf} -Oz -o ${label}.hard_filter.vcf.gz
   tabix ${label}.hard_filter.vcf.gz
   """ 

}
process Filter {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "3h"

   input:
   tuple val(label),file(vcf), file(vcf_index) from vcfs2

   output:
   tuple val(label),file("${label}.hard_filter.vcf.gz"), file("${label}.hard_filter.vcf.gz.tbi") into hard_filtered_vcfs2

   publishDir "${params.resultFolder}/vcfs", pattern: "${label}.hard_filter.vcf.gz*", mode: "move"

   """
   mkdir -p ${params.resultFolder}/vcfs
   bcftools filter -s FAIL -i '(TYPE="snp" & INFO/QD>=2 & INFO/FS<=60 & INFO/MQ>=40 & INFO/SOR<=3 & (INFO/MQRankSum="." | INFO/MQRankSum>=-12.4) & (INFO/ReadPosRankSum="." | INFO/ReadPosRankSum>=-8.0)) | (TYPE="indel" & INFO/QD>=2 & INFO/FS<=200 & INFO/SOR<=10 & (INFO/ReadPosRankSum="." | INFO/ReadPosRankSum>=-20))' ${vcf} -Oz -o ${label}.hard_filter.vcf.gz
   tabix ${label}.hard_filter.vcf.gz
   """ 

}

