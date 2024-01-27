#!/usr/bin/env nextflow

/*
* AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>; Daniel Taliun, PhD <daniel.taliun@mcgill.ca>
* VERSION: 2.0
* YEAR: 2023
*/

process GenotypeGVCFs {
   errorStrategy 'retry'
   maxRetries 0
   cache "lenient"
   cpus 1
   scratch '$SLURM_TMPDIR'
   stageInMode "copy"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref"
  
   input:
   tuple val(region), path(genomicsdb)

   output:
   tuple val(region), path("${region}.vcf.gz"), path("${region}.vcf.gz.tbi")

   """
   path=${region.split("_")[1]}_${region.split("_")[2]}
   gatk --java-options "-Xms4g -Xmx12g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs -OVI true -R /ref/${params.referenceGenome} -V gendb://\${path} -O ${region}.temp.vcf.gz --verbosity ERROR
   gatk --java-options "-Xmx8g -Xms4g" VariantFiltration -V ${region}.temp.vcf.gz -OVI true --filter-expression "ExcessHet > 54.69" --filter-name ExcessHet -O ${region}.vcf.gz --verbosity ERROR
   """
}


process GenotypeGVCFsWithPedigree {
   errorStrategy 'retry'
   maxRetries 0
   cache "lenient"
   cpus 1
   scratch '$SLURM_TMPDIR'
   stageInMode "copy"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref"
  
   input:
   tuple val(region), path(genomicsdb)
   each path(pedigree)

   output:
   tuple val(region), path("${region}.vcf.gz"), path("${region}.vcf.gz.tbi")

   """
   path=${region.split("_")[1]}_${region.split("_")[2]}
   gatk --java-options "-Xms2g -Xmx4g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs -OVI true -R /ref/${params.referenceGenome} --pedigree ${pedigree} -V gendb://\${path} -O ${region}.vcf.gz --verbosity ERROR
   """
}


process MakeSitesOnlyVcf {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1 
   scratch '$SLURM_TMPDIR'
   stageInMode "copy"

   container "${params.gatkContainer}"

   input:
   tuple val(region), path(vcf), path(vcf_index)
   
   output:
   tuple val(region), path("${region}.sitesonly.vcf.gz"), path("${region}.sitesonly.vcf.gz.tbi")

   """
   gatk --java-options "-Xms2g -Xmx4g" MakeSitesOnlyVcf -I ${vcf} -O ${region}.sitesonly.vcf.gz --VERBOSITY ERROR
   tabix -f -p vcf ${region}.sitesonly.vcf.gz
   """
}

process MergeSitesOnlyVCFs{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   scratch '$SLURM_TMPDIR'
   stageInMode "copy"

   input:   
   path vcfs_and_indices

   output:
   tuple path("all.sitesonly.vcf.gz"), path("all.sitesonly.vcf.gz.tbi")

   """
   find . -name "*.vcf.gz" | sort > files.txt
   bcftools concat -f files.txt -a -D -Oz -o all.sitesonly.vcf.gz
   tabix -f -p vcf all.sitesonly.vcf.gz
   """
}

process RecalibrateIndels {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   scratch '$SLURM_TMPDIR'
   stageInMode "copy"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.gatkBundle}:/bundle"

   input:
   tuple path(vcf_sites_only), path(vcf_sites_only_index)

   output:
   tuple path("indels.recal"), path("indels.recal.idx"), path("indels.tranches"), path("recalibrate_indels.output.plots.R")

   """
   gatk --java-options "-Xmx24g -Xms12g" VariantRecalibrator \
      -V ${vcf_sites_only} \
      --trust-all-polymorphic \
      -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
      -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
      -mode INDEL \
      --max-gaussians 4 \
      --resource:mills,known=false,training=true,truth=true,prior=12.0 /bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
      --resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 /bundle/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /bundle/Homo_sapiens_assembly38.dbsnp138.vcf \
      -O indels.recal \
      --tranches-file indels.tranches \
      --rscript-file recalibrate_indels.output.plots.R
   gatk IndexFeatureFile -I indels.recal 
   """
}


process RecalibrateSNVs {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   scratch '$SLURM_TMPDIR'
   stageInMode "copy"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.gatkBundle}:/bundle"

   input:
   tuple path(vcf_sites_only), path(vcf_sites_only_index)

   output:
   tuple path("snvs.recal"), path("snvs.recal.idx"), path("snvs.tranches"), path("recalibrate_snvs.output.plots.R")

   """
   gatk --java-options "-Xmx24g -Xms12g" VariantRecalibrator \
      -V ${vcf_sites_only} \
      --trust-all-polymorphic \
      -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
      -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP \
      -mode SNP \
      --max-gaussians 6 \
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /bundle/hapmap_3.3.hg38.vcf.gz \
      --resource:omni,known=false,training=true,truth=true,prior=12.0 /bundle/1000G_omni2.5.hg38.vcf.gz \
      --resource:1000G,known=false,training=true,truth=false,prior=10.0 /bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /bundle/Homo_sapiens_assembly38.dbsnp138.vcf \
      -O snvs.recal \
      --tranches-file snvs.tranches \
      --rscript-file recalibrate_snvs.output.plots.R
   gatk IndexFeatureFile -I snvs.recal
   """
}


process ApplyVQSR {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   scratch '$SLURM_TMPDIR'
   stageInMode "copy"

   container "${params.gatkContainer}"
  
   input:
   tuple val(region), path(vcf), path(vcf_index)
   tuple path(recal_table_snvs), path(recal_table_idx_snvs), path(recal_tranches_snvs), path(recal_plots_snvs)
   tuple path(recal_table_indels), path(recal_table_idx_indels), path(recal_tranches_indels), path(recal_plots_indels)

   output:
   tuple path("${region}.vqsr.vcf.gz"), path("${region}.vqsr.vcf.gz.tbi")

   storeDir "VCFs/${region.split('_')[0]}/"
   
   """
   gatk --java-options "-Xmx12g -Xms4g" ApplyVQSR -OVI true -V ${vcf} --recal-file ${recal_table_indels} --tranches-file ${recal_tranches_indels} --truth-sensitivity-filter-level 99.9 -mode INDEL -O temp.vcf.gz
   gatk --java-options "-Xmx12g -Xms4g" ApplyVQSR -OVI true -V temp.vcf.gz --recal-file ${recal_table_snvs} --tranches-file ${recal_tranches_snvs} --truth-sensitivity-filter-level 99.9 -mode SNP -O ${region}.vqsr.vcf.gz
   """
}


workflow {
	chunked_genomicsDB = Channel.fromPath(params.genomicsDB, type: "dir").map{ file-> [ file.getParent().toString().split('/').last() + "_" + file.getName(), file] }
	if (params.pedigree) {
		chunked_vcfs = GenotypeGVCFsWithPedigree(chunked_genomicsDB, Channel.fromPath(params.pedigree))
	} else {
		chunked_vcfs = GenotypeGVCFs(chunked_genomicsDB)
	}
	sites_only_chunked_vcfs = MakeSitesOnlyVcf(chunked_vcfs)
	sites_only_vcfs = MergeSitesOnlyVCFs(sites_only_chunked_vcfs.map{ [it[1], it[2]] }.collect())
	recal_table_tranches_indels = RecalibrateIndels(sites_only_vcfs)
	recal_table_tranches_snvs = RecalibrateSNVs(sites_only_vcfs)
	ApplyVQSR(chunked_vcfs, recal_table_tranches_snvs, recal_table_tranches_indels)
}
