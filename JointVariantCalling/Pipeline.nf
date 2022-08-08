#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2021
*/

inputFiles = Channel.fromPath(params.inputFiles).map { file -> [file.getBaseName(), file] }
intervals = Channel.fromPath(params.intervals).map { file -> [file.getBaseName(), file]}
regions = Channel.from(file(params.intervalsList).readLines())
chroms = Channel.from(file(params.chroms).readLines())

process HaplotypeCaller {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref"

   input:   
   tuple val(interval_label), file(interval_list), val(input_label), file(input) from intervals.combine(inputFiles)

   output:
   tuple val(input_label), file("${interval_label}.${input_label}.g.vcf.gz"), file("${interval_label}.${input_label}.g.vcf.gz.tbi") into interval_vcfs

   when:
   params.runHaplotypeCaller
   
   """
   gatk --java-options -Xmx14G HaplotypeCaller --ERC GVCF --native-pair-hmm-threads 1 -G StandardAnnotation -G StandardHCAnnotation -L ${interval_list} -R /ref/${params.referenceGenome}.fa -I ${input} -O ${interval_label}.${input_label}.g.vcf.gz
   """
}


process Merge {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "12h"

   input:
   tuple val(label), file(vcfs), file(vcf_indices) from interval_vcfs.groupTuple(by: 0)

   output:
   file("${label}.g.vcf.gz*") into gvcfs

   publishDir "${params.gvcfsResultFolder}", pattern: "${label}.g.vcf*", mode: "move"

   """
   find . -name "*.g.vcf.gz" | sort > files.txt
   bcftools concat -a -f files.txt -Oz -o ${label}.g.vcf.gz
   tabix ${label}.g.vcf.gz
   """
}


process GenomicsDBImport {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB" 
   time "48h"
   scratch '$SLURM_TMPDIR'

   container "${params.gatkContainer}"
   containerOptions "-B ${params.genomicsDBImportFolder}:/GenomicsDBImportFolder -B ${params.genomicsDB}:/GenomicsDB -B ${params.projectFolder}:/samples -B ${params.gvcfsResultFolder}:/input -B ${params.referenceDir}:/ref"

   input:   
   val(chr) from regions

   when:
   params.runGenomicDBImport

   script:
   if ( params.genomicDBImportMode == "Initialize" )
      """ 
      mkdir -p region/${chr}
      gatk --java-options "-Xms4g -Xmx12g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport --genomicsdb-workspace-path region/${chr}/GenomicsDBImport --batch-size 200 --intervals ${chr} --sample-name-map /samples/BQC19SampleMap.txt
      mkdir -p /GenomicsDBImportFolder/${chr}
      cp -R region/${chr}/* /GenomicsDBImportFolder/${chr}
      """
   else if ( params.genomicDBImportMode == "Add" )
      """ 
      mkdir -p region
      cp -r /GenomicsDB/${chr} region/
      gatk --java-options "-Xms4g -Xmx12g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport --genomicsdb-update-workspace-path region/${chr}/GenomicsDBImport --batch-size 200 --intervals ${chr} --sample-name-map /samples/BQC19SampleMap_5.txt
      mkdir -p /GenomicsDBImportFolder/${chr}
      cp -R region/${chr}/* /GenomicsDBImportFolder/${chr}
      """
   

}

process GenotypeGVCFs {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB" 
   time "96h"
   scratch '$SLURM_TMPDIR'

   container "${params.gatkContainer}"
   containerOptions "-B ${params.genomicsDBImportFolder}:/GenomicsDBImportFolder  -B ${params.jvcMSResultFolder}:/vcfs -B ${params.referenceDir}:/ref"

   input:   
   val(chr) from regions

   when:
   params.runGenomicDBImport
   
   """
   mkdir work
   cp -r /GenomicsDBImportFolder/${chr}/ ./work
   cp /ref/${params.referenceGenome}.fa ./work
   cp /ref/${params.referenceGenome}.fa.fai ./work
   cp /ref/${params.referenceGenome}.dict ./work
   gatk --java-options "-Xms4g -Xmx12g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs -R ./work/${params.referenceGenome}.fa -V gendb://work/${chr}/GenomicsDBImport/ -O ./work/${chr}.jvc.vcf.gz --verbosity ERROR
   cp ./work/${chr}.jvc.vcf.gz /vcfs/
   """

}

process VariantFilter{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB" 
   time "12h"
   scratch '$SLURM_TMPDIR'

   container "${params.gatkContainer}"
   containerOptions "-B ${params.jvcMSResultFolder}:/vcfs -B ${params.vqsrFolder}:/vqsr"

   input:   
   val(vcfFile) from regions
   
   output:
   val "${vcfFile}" into excesshetFiltered

   when:
   params.runVariantFilter

   """
   cp /vcfs/${vcfFile}.jvc.vcf.gz ./
   tabix -f -p vcf ./${vcfFile}.jvc.vcf.gz
   gatk --java-options "-Xmx8g -Xms4g" VariantFiltration -V ./${vcfFile}.jvc.vcf.gz -OVI true --filter-expression "ExcessHet > 54.69" --filter-name ExcessHet -O ./${vcfFile}.excesshet.jvc.vcf.gz --verbosity ERROR
   cp ./${vcfFile}.excesshet.* /vqsr
   """
}

process MakeSitesOnlyVcf{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB" 
   time "6h"
   scratch '$SLURM_TMPDIR'

   container "${params.gatkContainer}"
   containerOptions "-B ${params.vqsrFolder}:/vqsr"

   input:   
   val vcfFile from excesshetFiltered
   
   output:
   val "${vcfFile}" into sitesOnlyVcf

   """
   cp /vqsr/${vcfFile}.excesshet* ./
   gatk MakeSitesOnlyVcf -I ${vcfFile}.excesshet.jvc.vcf.gz -O ${vcfFile}.sitesonly.excesshet.jvc.vcf.gz --VERBOSITY ERROR
   tabix -f -p vcf ${vcfFile}.sitesonly.excesshet.jvc.vcf.gz
   cp ./${vcfFile}.sitesonly.excesshet.jvc.vcf.gz* /vqsr
   """
}

process MergeSiteonlyVCFs{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "64 GB" 
   time "6h"

   input:   
   val(vcf) from sitesOnlyVcf.collect()

   output:
   val true into variantrecal

   """
   find ${params.vqsrFolder} -name "*.sitesonly.excesshet.jvc.vcf.gz" | sort > files.txt
   bcftools concat -f files.txt -a -D -Oz -o ${params.vqsrFolder}/${params.label}.sitesonly.vcf.gz
   tabix -f -p vcf ${params.vqsrFolder}/${params.label}.sitesonly.vcf.gz
   """
}
variantrecal.into { vr_Indels; vr_SNPs }
process VariantRecalibratorIndels{
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "64 GB" 
   time "24h"
   scratch '$SLURM_TMPDIR'

   container "${params.gatkContainer}"
   containerOptions "-B ${params.vqsrFolder}:/vqsr -B ${params.gatkBundle}:/bundle"

   input:   
   val flag from vr_Indels
   
   output:
   val true into vrIndel

   """
   cp /vqsr/BQC19.sitesonly.vcf.gz* ./
   cp /bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf* ./
   cp /bundle/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf* ./
   cp /bundle/Homo_sapiens_assembly38.dbsnp138* ./
   gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator -V ./${params.label}.sitesonly.vcf.gz --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -mode INDEL --max-gaussians 4 --resource:mills,known=false,training=true,truth=true,prior=12.0 ./Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 ./Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./Homo_sapiens_assembly38.dbsnp138.vcf -O ./${params.label}_indels.recal --tranches-file ./${params.label}_indels.tranches --rscript-file ./IndelRecal.output.plots.R
   gatk IndexFeatureFile -I ./${params.label}_indels.recal
   cp ./${params.label}_indels.recal* /vqsr
   cp ./${params.label}_indels.tranches /vqsr
   cp ./IndelRecal.output.plots.R /vqsr
   """
}

process VariantRecalibratorSNPs{
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "64 GB" 
   time "24h"
   scratch '$SLURM_TMPDIR'

   container "${params.gatkContainer}"
   containerOptions "-B ${params.vqsrFolder}:/vqsr -B ${params.gatkBundle}:/bundle"

   input:   
   val flag from vr_SNPs
   
   output:
   val true into vrSNP

   """
   cp /vqsr/BQC19.sitesonly.vcf.gz* ./
   cp /bundle/hapmap_3.3.hg38.vcf.gz* ./
   cp /bundle/1000G_omni2.5.hg38.vcf.gz* ./
   cp /bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz* ./
   cp /bundle/Homo_sapiens_assembly38.dbsnp138* ./
   gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator -V ./${params.label}.sitesonly.vcf.gz --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP -mode SNP --max-gaussians 6 --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ./hapmap_3.3.hg38.vcf.gz --resource:omni,known=false,training=true,truth=true,prior=12.0 ./1000G_omni2.5.hg38.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 ./1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ./Homo_sapiens_assembly38.dbsnp138.vcf -O ./${params.label}_snps.recal --tranches-file ./${params.label}_snps.tranches --rscript-file ./SNPRecal.output.plots.R
   gatk IndexFeatureFile -I ./${params.label}_snps.recal
   cp ./${params.label}_snps.recal* /vqsr
   cp ./${params.label}_snps.tranches /vqsr
   cp ./SNPRecal.output.plots.R /vqsr
   """
}

process ApplyVQSR{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB" 
   time "24h"
   scratch '$SLURM_TMPDIR'

   container "${params.gatkContainer}"
   containerOptions "-B ${params.vqsrFolder}:/vqsr"

   input:   
   val(vcfFile) from regions
   val(indels) from vrIndel.collect()
   val(snps) from vrSNP.collect()


   output:
   val "${vcfFile}" into applyVQSRSNPVCFs
   
   """
   cp /vqsr/${vcfFile}.excesshet.jvc.vcf.gz* ./
   cp /vqsr/${params.label}_indels.recal* ./
   cp /vqsr/${params.label}_indels.tranches ./
   cp /vqsr/${params.label}_snps.recal* ./
   cp /vqsr/${params.label}_snps.tranches ./
   gatk --java-options "-Xmx12g -Xms4g" ApplyVQSR -V ./${vcfFile}.excesshet.jvc.vcf.gz --recal-file ./${params.label}_indels.recal --tranches-file ./${params.label}_indels.tranches --truth-sensitivity-filter-level 99.9 --create-output-variant-index true -mode INDEL -O ./${vcfFile}.indel.recalibrated.vcf.gz
   gatk --java-options "-Xmx12g -Xms4g" ApplyVQSR -V ./${vcfFile}.indel.recalibrated.vcf.gz --recal-file ./${params.label}_snps.recal --tranches-file ./${params.label}_snps.tranches --truth-sensitivity-filter-level 99.9 --create-output-variant-index true -mode SNP -O ./${vcfFile}.vqsr.vcf.gz
   cp ./${vcfFile}.vqsr.vcf.gz* /vqsr
   """
}


process MergeVCFsByChrom{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB" 
   time "2h"

   input:
   val(chr) from chroms
   val(vcfs) from applyVQSRSNPVCFs.collect()

   output:
   tuple val(chr), file("${chr}_files.txt") into mergeFiles

   """
   find ${params.vqsrFolder} -name "${chr}:*.vqsr.vcf.gz" | sort > ${chr}_files.txt
   """

}

process MergeFiles{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB" 
   time "24h"
  
   container "${params.gatkContainer}"
   containerOptions "-B ${params.result}:/result"

   input:   
   tuple val(chr), file(files) from mergeFiles

   output:
   file("${params.label}.${chr}.sitesonly.vqsr.vcf.gz*") into sitesOnlyVQSRVcf

   """
   bcftools concat -f ${files} -a -D -Oz -o /result/${params.label}.${chr}.vqsr.vcf.gz
   tabix -f -p vcf /result/${params.label}.${chr}.vqsr.vcf.gz
   gatk MakeSitesOnlyVcf -I  /result/${params.label}.${chr}.vqsr.vcf.gz -O ${params.label}.${chr}.sitesonly.vqsr.vcf.gz --VERBOSITY ERROR
   tabix -f -p vcf ${params.label}.${chr}.sitesonly.vqsr.vcf.gz
   """
}

process MergeSiteonlyVQSRVCFs{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB" 
   time "24h"

   input:   
   file('*') from sitesOnlyVQSRVcf.collect()

   """
   find . -name "*.sitesonly.vqsr.vcf.gz" | sort > files.txt
   bcftools concat -f files.txt -a -D -Oz -o ${params.result}/${params.label}.sitesonly.vqsr.vcf.gz
   tabix -f -p vcf ${params.result}/${params.label}.sitesonly.vqsr.vcf.gz
   """
}