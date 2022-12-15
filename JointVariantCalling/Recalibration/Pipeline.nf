#!/usr/bin/env nextflow

/*
* AUTHORS: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>; Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2022
*/


process BaseRecalibrator {
   errorStrategy 'finnish'
   cache "lenient"
   cpus 1
   memory "5 GB"
   time "12h"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref -B ${params.bundle}:/bundle"

   input:
   tuple path(input), val(index_type)
   
   output:
   tuple path(input), val(index_type), path("${input.getBaseName()}.table")

   script:
   if (( params.bundleBuild == "GRCh37" ) || ( params.bundleBuild == "hg19" ))
      """
      gatk --java-options "-Xmx4G" BaseRecalibrator -R /ref/${params.referenceGenome} --known-sites /bundle/Homo_sapiens_assembly19.dbsnp138.vcf --known-sites /bundle/Homo_sapiens_assembly19.known_indels.vcf --known-sites /bundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz -I ${input} -O ${input.getBaseName()}.table
      """
   else if (( params.bundleBuild == "GRCh38") || ( params.bundleBuild == "hg38" ))
      """
      gatk --java-options "-Xmx4G" BaseRecalibrator -R /ref/${params.referenceGenome} --known-sites /bundle/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /bundle/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -I ${input} -O ${input.getBaseName()}.table
      """
   else
   	error "Invalid GATK bundle assembly name: ${params.bundleBuild}"
}


process ApplyBQSR {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "5 GB"
   time "12h"
   scratch '$SLURM_TMPDIR'

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref"

   input:
   tuple path(input), val(index_type), path(recal_table)

   output:
   tuple path("${input.getBaseName()}.recal.${input.getExtension()}"), path("${input.getBaseName()}.recal.${index_type}")

   publishDir "final", pattern: "${input.getBaseName()}.recal.*"

   """
   gatk --java-options "-Xmx4G" ApplyBQSR --create-output-bam-index false -I ${input} -R /ref/${params.referenceGenome} --bqsr-recal-file ${recal_table} -O ${input.getBaseName()}.recal.${input.getExtension()}
   samtools index ${input.getBaseName()}.recal.${input.getExtension()}
   """
}


workflow {
	bams_and_recal_table = BaseRecalibrator(Channel.fromPath(params.inputFiles).map { file -> [file, file.getExtension() == "bam" ? "bai" : "crai"] })
	ApplyBQSR(bams_and_recal_table)	
}
