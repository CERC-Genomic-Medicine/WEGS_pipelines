#!/usr/bin/env nextflow

/*
* AUTHORS: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>; Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2021
*/


process BaseRecalibrator {
   errorStrategy "finish"
   cache "lenient"
   cpus 1
   memory "5 GB"
   time "12h"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref -B ${params.bundle}:/bundle"

   input:   
   tuple path(input), val(index_type)

   output:
   tuple path("${input.getBaseName()}.recal.${input.getExtension()}"), path("${input.getBaseName()}.recal.${index_type}")

   publishDir "${params.resultFolder}/bams", pattern: "${input.getBaseName()}.recal.*", mode: "copy"

   script:
   if (( params.bundleBuild == "GRCh37" ) || ( params.bundleBuild == "hg19" ))
   	"""
   	gatk --java-options "-Xmx4G" BaseRecalibrator -R /ref/${params.referenceGenome} --known-sites /bundle/Homo_sapiens_assembly19.dbsnp138.vcf --known-sites /bundle/Homo_sapiens_assembly19.known_indels.vcf --known-sites /bundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz -I ${input} -O ${input.getBaseName()}.table
	gatk --java-options "-Xmx4G" ApplyBQSR --create-output-bam-index true -I ${input} -R /ref/${params.referenceGenome} --bqsr-recal-file ${input.getBaseName()}.table -O ${input.getBaseName()}.recal.${input.getExtension()}
	"""
   else if (( params.bundleBuild == "GRCh38") || ( params.bundleBuild == "hg38" ))
   	"""
   	gatk --java-options "-Xmx4G" BaseRecalibrator -R /ref/${params.referenceGenome} --known-sites /bundle/Homo_sapiens_assembly38.dbsnp138.vcf --known-sites /bundle/Homo_sapiens_assembly38.known_indels.vcf.gz --known-sites /bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz -I ${input} -O ${input.getBaseName()}.table
	gatk --java-options "-Xmx4G" ApplyBQSR --create-output-bam-index true -I ${input} -R /ref/${params.referenceGenome} --bqsr-recal-file ${input.getBaseName()}.table -O ${input.getBaseName()}.recal.${input.getExtension()}
   	"""
   else
   	error "Invalid GATK bundle assembly name: ${params.bundleBuild}"
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
   tuple val(interval_label), path(interval_list), path(input), path(index)

   output:
   tuple path("${interval_label}.${input.getSimpleName()}.vcf.gz"), path("${interval_label}.${input.getSimpleName()}.vcf.gz.tbi")

   """
   gatk --java-options -Xmx12G HaplotypeCaller --native-pair-hmm-threads 1 -G StandardAnnotation -G StandardHCAnnotation -L ${interval_list} -R /ref/${params.referenceGenome} -I ${input} -O ${interval_label}.${input.getSimpleName()}.vcf.gz
   """
}


process MergeAndNormalize {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "6h"

   input:
   tuple val(label), path(vcfs), path(vcf_indices)

   output:
   tuple(path("${label}.norm.vcf.gz"), path("${label}.norm.vcf.gz.tbi"), emit: merged_normalized_vcfs)
   tuple(path("${label}.vcf.gz"), path("${label}.vcf.gz.tbi"), emit: merged_vcfs)

   publishDir "${params.resultFolder}/VCF_norm_nofilter", pattern: "${label}.norm.vcf.gz*", mode: "copy"
   publishDir "${params.resultFolder}/VCF_raw", pattern: "${label}.vcf.gz*", mode: "copy"

   """
   find . -name "*.vcf.gz" | sort > files.txt
   bcftools concat -f files.txt -Oz -o ${label}.vcf.gz
   tabix ${label}.vcf.gz

   bcftools norm -f ${params.referenceDir}/${params.referenceGenome} -m - ${label}.vcf.gz | bcftools view -v snps,indels -Oz -o ${label}.norm.vcf.gz
   tabix ${label}.norm.vcf.gz
   """
}


process HardFilter {
   errorStrategy 'finish'
   maxRetries 1
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "3h"

   container "${params.gatkContainer}"

   input:
   tuple val(name), path(vcf), path(vcf_index)

   output:
   tuple path("${name}.hard_filter.vcf.gz"), path("${name}.hard_filter.vcf.gz.tbi")

   publishDir "${params.resultFolder}/VCF_norm_hard_filter", pattern: "${name}.hard_filter.vcf.gz*", mode: "copy"

   """
   gatk --java-options -Xmx4G SelectVariants -V ${vcf} -select-type SNP -O snps.vcf.gz
   gatk --java-options -Xmx4G SelectVariants -V ${vcf} -select-type INDEL -O indels.vcf.gz


   gatk --java-options -Xmx4G VariantFiltration -V snps.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "SOR > 3.0" --filter-name "SOR3" -filter "FS > 60.0" --filter-name "FS60" -filter "MQ < 40.0" --filter-name "MQ40" -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" -O snps.filtered.vcf.gz

   gatk --java-options -Xmx4G VariantFiltration -V indels.vcf.gz -filter "QD < 2.0" --filter-name "QD2" -filter "QUAL < 30.0" --filter-name "QUAL30" -filter "FS > 200.0" --filter-name "FS200" -filter "ReadPosRankSum < -20.0" --filter-name "ReadPosRankSum-20" -O indels.filtered.vcf.gz

   gatk --java-options -Xmx4G MergeVcfs -I snps.filtered.vcf.gz -I indels.filtered.vcf.gz -O ${name}.hard_filter.vcf.gz
   """ 
}


process CNNFilter1D {
   errorStrategy 'finish'
   maxRetries 1
   cache "lenient"
   cpus 1
   memory '8 GB'
   time '24h'

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref -B ${params.bundle}:/bundle"

   input:
   tuple val(name), path(vcf), path(vcf_index)

   output:
   tuple path("${name}.cnn1d_filter.vcf.gz"), path("${name}.cnn1d_filter.vcf.gz.tbi")

   publishDir "${params.resultFolder}/VCF_norm_cnn1d_filter", pattern: "${name}.cnn1d_filter.vcf.gz*", mode: "copy"

   script:
   if (( params.bundleBuild == "GRCh37" ) || ( params.bundleBuild == "hg19" ))
   """
      gatk CNNScoreVariants -R /ref/${params.referenceGenome} -V ${vcf} -O annotated.vcf.gz
      gatk FilterVariantTranches --resource /bundle/1000G_omni2.5.b37.vcf.gz --resource /bundle/hapmap_3.3.b37.vcf.gz --resource /bundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz --info-key CNN_1D --snp-tranche 99.95 --indel-tranche 99.4 -V annotated.vcf.gz -O ${name}.cnn1d_filter.vcf.gz
   """
   else if (( params.bundleBuild == "GRCh38") || ( params.bundleBuild == "hg38" ))
   """
      gatk CNNScoreVariants -R /ref/${params.referenceGenome} -V ${vcf} -O annotated.vcf.gz
      gatk FilterVariantTranches --resource /bundle/1000G_omni2.5.hg38.vcf.gz --resource /bundle/hapmap_3.3.hg38.vcf.gz --resource /bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --info-key CNN_1D --snp-tranche 99.95 --indel-tranche 99.4 -V annotated.vcf.gz -O ${name}.cnn1d_filter.vcf.gz
   """
   else
      error "Invalid GATK bundle assembly name: ${params.bundleBuild}"
}


process CNNFilter2D {
   errorStrategy "finish"
   maxRetries 1
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "24h"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref -B ${params.bundle}:/bundle"

   input:
   tuple val(name), path(vcf), path(vcf_index), path(bam), path(bam_index)

   output:
   tuple path("${name}.cnn2d_filter.vcf.gz"), path("${name}.cnn2d_filter.vcf.gz.tbi")

   publishDir "${params.resultFolder}/VCF_norm_cnn2d_filter", pattern: "${name}.cnn2d_filter.vcf.gz*", mode: "copy"

   script:
   if (( params.bundleBuild == "GRCh37" ) || ( params.bundleBuild == "hg19" ))
      """
      gatk CNNScoreVariants --tensor-type read_tensor -R /ref/${params.referenceGenome} -I ${bam} -V ${vcf} -O annotated.vcf.gz
      gatk FilterVariantTranches --resource /bundle/1000G_omni2.5.b37.vcf.gz --resource /bundle/hapmap_3.3.b37.vcf.gz --resource /bundle/Mills_and_1000G_gold_standard.indels.b37.vcf.gz --info-key CNN_2D --snp-tranche 99.95 --indel-tranche 99.4 -V annotated.vcf.gz -O ${name}.cnn2d_filter.vcf.gz
      """
   else if (( params.bundleBuild == "GRCh38") || ( params.bundleBuild == "hg38" ))
      """
      gatk CNNScoreVariants --tensor-type read_tensor -R /ref/${params.referenceGenome} -I ${bam} -V ${vcf} -O annotated.vcf.gz
      gatk FilterVariantTranches --resource /bundle/1000G_omni2.5.hg38.vcf.gz --resource /bundle/hapmap_3.3.hg38.vcf.gz --resource /bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --info-key CNN_2D --snp-tranche 99.95 --indel-tranche 99.4 -V annotated.vcf.gz -O ${name}.cnn2d_filter.vcf.gz
      """
   else
      error "Invalid GATK bundle assembly name: ${params.bundleBuild}"
}


workflow {
	intervals = Channel.fromPath(params.intervals).map { file -> [file.getParent().toString().split('/').last(), file] }

	if (params.doBQSR == true) { // perform base recallibration if needed
		BaseRecalibrator(Channel.fromPath(params.inputFiles).map { file -> [file, file.getExtension() == "bam" ? "bai" : "crai"] })
		bams = BaseRecalibrator.out
	} else {
		bams = Channel.fromPath(params.inputFiles).map { file -> [file, file + (file.getExtension() == "bam" ? ".bai" : ".crai")] }
	}

	chunked_vcfs = HaplotypeCaller(intervals.combine(bams))

	MergeAndNormalize(chunked_vcfs.map { it -> [ it[0].getName().toString().split('\\.')[1], it[0], it[1] ]  }.groupTuple(by: 0))
	merged_normalized_vcfs = MergeAndNormalize.out.merged_normalized_vcfs.map { it -> [ it[0].getName().toString().replace(".vcf.gz", ""), it[0], it[1] ] }
	HardFilter(merged_normalized_vcfs)
	CNNFilter1D(merged_normalized_vcfs)

	bams_with_key = bams.map { list -> [list[0].getSimpleName(), list[0], list[1]] }
	merged_normalized_vcfs_with_key = merged_normalized_vcfs.map { list -> [list[1].getSimpleName(), list[0], list[1], list[2] ] }
	merged_normalized_vcfs_and_bams = merged_normalized_vcfs_with_key.join(bams_with_key, by: 0).map { list -> [list[1], list[2], list[3], list[4], list[5]] }
	CNNFilter2D(merged_normalized_vcfs_and_bams)
}
