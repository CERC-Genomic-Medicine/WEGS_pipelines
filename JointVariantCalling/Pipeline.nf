#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2021
*/

inputFiles = Channel.fromPath(params.inputFiles).map { file -> [file.getBaseName(), file] }
intervals = Channel.fromPath(params.intervals).map { file -> [file.getBaseName(), file]}
regions = Channel.from(file(params.intervalsList).readLines())

process HaplotypeCaller {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "5 GB"
   time "8h"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.referenceDir}:/ref"

   input:   
   tuple val(interval_label), file(interval_list), val(input_label), file(input) from intervals.combine(inputFiles)

   output:
   tuple val(input_label), file("${interval_label}.${input_label}.g.vcf.gz"), file("${interval_label}.${input_label}.g.vcf.gz.tbi") into interval_vcfs

   """
   gatk --java-options -Xmx4G HaplotypeCaller --ERC GVCF --native-pair-hmm-threads 1 -G StandardAnnotation -G StandardHCAnnotation -L ${interval_list} -R /ref/${params.referenceGenome}.fa -I ${input} -O ${interval_label}.${input_label}.g.vcf.gz
   """
}


process Merge {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "4 GB"
   time "8h"

   input:
   tuple val(label), file(vcfs), file(vcf_indices) from interval_vcfs.groupTuple(by: 0)

   output:
   file("${label}.g.vcf.gz*") into gvcfs

   publishDir "${params.gvcfsResultFolder}", pattern: "${label}.g.vcf*"

   """
   find . -name "*.g.vcf.gz" | sort > files.txt
   bcftools concat -a -f files.txt -Oz -o ${label}.g.vcf.gz
   tabix ${label}.g.vcf.gz
   """
}


process GenotypeGVCFs {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "64 GB"
   time "48h"
   maxForks 1

   container "${params.gatkContainer}"
   containerOptions "-B ${params.tempFolder}:/tmp -B ${params.genomicsDBImportFolder}:/GenomicsDBImportFolder -B ${params.sampleMapFolder}:/samples -B ${params.gvcfsResultFolder}:/input -B ${params.referenceDir}:/ref"

   input:   
   file('*') from gvcfs.collect()

   output:
   file("output.vcf.gz") into jvcVcf
   
   publishDir "${params.jvcMSResultFolder}"
   
   """
   gatk --java-options "-Xms4g -Xmx48g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport --genomicsdb-workspace-path /GenomicsDBImportFolder/GenomicsDBImport --batch-size 50 --tmp-dir /tmp --intervals /ref/${params.referenceGenome}.intervals.bed --sample-name-map /samples/BQC19SampleMap.txt
   gatk --java-options "-Xms4g -Xmx48g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs -R /ref/${params.referenceGenome}.fa -V gendb:///GenomicsDBImportFolder/GenomicsDBImport/ -O BQC19.jvc.vcf.gz
   """

}

/*process GenomicsDBImport {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "64 GB" 
   time "12h"
   scratch true

   container "${params.gatkContainer}"
   containerOptions "-B ${params.genomicsDBImportFolder}:/GenomicsDBImportFolder -B ${params.sampleMapFolder}:/samples -B ${params.gvcfsResultFolder}:/input -B ${params.referenceDir}:/ref"

   input:   
   val(chr) from regions
   
   """
   if [ -d "/GenomicsDBImportFolder/regionWise/${chr}" ]; then rm -Rf /GenomicsDBImportFolder/regionWise/${chr}; fi
   mkdir $TMPDIR/${chr}
   gatk --java-options "-Xms4g -Xmx60g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport --genomicsdb-workspace-path $TMPDIR/${chr}/GenomicsDBImport --batch-size 200 --tmp-dir $TMPDIR --intervals ${chr} --sample-name-map /samples/BQC19SampleMap.txt
   mkdir /GenomicsDBImportFolder/regionWise/${chr}
   cp -R $TMPDIR/${chr}/* /GenomicsDBImportFolder/regionWise/${chr}
   """

}*/

/*process GenotypeGVCFs {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "80 GB" 
   time "48h"
   scratch true

   container "${params.gatkContainer}"
   containerOptions "-B ${params.genomicsDBImportFolder}:/GenomicsDBImportFolder  -B ${params.jvcMSResultFolder}:/vcfs -B ${params.referenceDir}:/ref"

   input:   
   val(chr) from regions
   
   """
   cp -r /GenomicsDBImportFolder/regionWise/${chr}/ $TMPDIR/
   cp /ref/${params.referenceGenome}.fa $TMPDIR/
   cp /ref/${params.referenceGenome}.fa.fai $TMPDIR/
   cp /ref/${params.referenceGenome}.dict $TMPDIR/
   gatk --java-options "-Xms4g -Xmx32g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenotypeGVCFs -R $TMPDIR/${params.referenceGenome}.fa -V gendb:///$TMPDIR/${chr}/GenomicsDBImport/ -O $TMPDIR/${chr}.jvc.vcf.gz --verbosity ERROR
   cp $TMPDIR/${chr}.jvc.vcf.gz /vcfs/
   """

}*/

//vcfFiles = Channel.fromList( ['chr1:10001-207666'] )
/*process VariantFiltration{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB" 
   time "4h"
   scratch true

   container "${params.gatkContainer}"
   containerOptions "-B ${params.jvcMSResultFolder}:/vcfs -B ${params.vqsrFolder}:/vqsr"

   input:   
   val(vcfFile) from regions
   
   output:
   val "${vcfFile}" into excesshetFiltered

   """
   cp /vcfs/${vcfFile}.jvc* $TMPDIR/
   gatk --java-options "-Xmx8g -Xms4g" VariantFiltration -V $TMPDIR/${vcfFile}.jvc.vcf.gz -OVI true --filter-expression "ExcessHet > 54.69" --filter-name ExcessHet -O file:////$TMPDIR/${vcfFile}.excesshet.jvc.vcf.gz --verbosity ERROR
   cp $TMPDIR/${vcfFile}.excesshet.* /vqsr
   """
}

process MakeSitesOnlyVcf{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB" 
   time "4h"
   scratch true

   container "${params.gatkContainer}"
   containerOptions "-B ${params.vqsrFolder}:/vqsr"

   input:   
   val vcfFile from excesshetFiltered
   
   output:
   val "${vcfFile}" into sitesOnlyVcf

   """
   cp /vqsr/${vcfFile}.excesshet* $TMPDIR/
   gatk MakeSitesOnlyVcf -I $TMPDIR/${vcfFile}.excesshet.jvc.vcf.gz -O $TMPDIR/${vcfFile}.sitesonly.excesshet.jvc.vcf.gz --VERBOSITY ERROR
   tabix -f -p vcf $TMPDIR/${vcfFile}.sitesonly.excesshet.jvc.vcf.gz
   cp $TMPDIR/${vcfFile}.sitesonly.* /vqsr
   """
}

process MergeSiteonlyVCFs{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "32 GB" 
   time "6h"

   input:   
   val(vcf) from sitesOnlyVcf.collect()

   output:
   val true into variantrecal

   """
   bcftools concat -f ${params.vqsrFolder}/list.txt -a -D -Oz -o ${params.vqsrFolder}/BQC19.sitesonly.vcf.gz
   tabix -f -p vcf ${params.vqsrFolder}/BQC19.sitesonly.vcf.gz
   """
}
variantrecal.into { vr_Indels; vr_SNPs }
process VariantRecalibratorIndels{
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "64 GB" 
   time "3h"
   scratch true

   container "${params.gatkContainer}"
   containerOptions "-B ${params.vqsrFolder}:/vqsr -B ${params.gatkBundle}:/bundle"

   input:   
   val flag from vr_Indels
   
   output:
   val true into vrIndel

   """
   cp /vqsr/BQC19.sitesonly.vcf.gz* $TMPDIR/
   cp /bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf* $TMPDIR/
   cp /bundle/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf* $TMPDIR/
   cp /bundle/Homo_sapiens_assembly38.dbsnp138* $TMPDIR/
   gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator -V $TMPDIR/BQC19.sitesonly.vcf.gz --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP -mode INDEL --max-gaussians 4 --resource:mills,known=false,training=true,truth=true,prior=12.0 $TMPDIR/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz --resource:axiomPoly,known=false,training=true,truth=false,prior=10.0 $TMPDIR/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $TMPDIR/Homo_sapiens_assembly38.dbsnp138.vcf -O $TMPDIR/BQC19_indels.recal --tranches-file $TMPDIR/BQC19_indels.tranches --rscript-file $TMPDIR/IndelRecal.output.plots.R
   gatk IndexFeatureFile -I $TMPDIR/BQC19_indels.recal
   cp $TMPDIR/BQC19_indels.recal* /vqsr
   cp $TMPDIR/BQC19_indels.tranches /vqsr
   cp $TMPDIR/IndelRecal.output.plots.R /vqsr
   """
}

process VariantRecalibratorSNPs{
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "64 GB" 
   time "3h"
   scratch true

   container "${params.gatkContainer}"
   containerOptions "-B ${params.vqsrFolder}:/vqsr -B ${params.gatkBundle}:/bundle"

   input:   
   val flag from vr_SNPs
   
   output:
   val true into vrSNP

   """
   cp /vqsr/BQC19.sitesonly.vcf.gz* $TMPDIR/
   cp /bundle/hapmap_3.3.hg38.vcf.gz* $TMPDIR/
   cp /bundle/1000G_omni2.5.hg38.vcf.gz* $TMPDIR/
   cp /bundle/1000G_phase1.snps.high_confidence.hg38.vcf.gz* $TMPDIR/
   cp /bundle/Homo_sapiens_assembly38.dbsnp138* $TMPDIR/
   gatk --java-options "-Xmx24g -Xms24g" VariantRecalibrator -V $TMPDIR/BQC19.sitesonly.vcf.gz --trust-all-polymorphic -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 -an QD -an MQRankSum -an ReadPosRankSum -an FS -an MQ -an SOR -an DP -mode SNP --max-gaussians 6 --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $TMPDIR/hapmap_3.3.hg38.vcf.gz --resource:omni,known=false,training=true,truth=true,prior=12.0 $TMPDIR/1000G_omni2.5.hg38.vcf.gz --resource:1000G,known=false,training=true,truth=false,prior=10.0 $TMPDIR/1000G_phase1.snps.high_confidence.hg38.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $TMPDIR/Homo_sapiens_assembly38.dbsnp138.vcf -O $TMPDIR/BQC19_snps.recal --tranches-file $TMPDIR/BQC19_snps.tranches --rscript-file $TMPDIR/SNPRecal.output.plots.R
   gatk IndexFeatureFile -I $TMPDIR/BQC19_snps.recal
   cp $TMPDIR/BQC19_snps.recal* /vqsr
   cp $TMPDIR/BQC19_snps.tranches /vqsr
   cp $TMPDIR/SNPRecal.output.plots.R /vqsr
   """
}

/*process ApplyVQSRIndel{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB" 
   time "4h"
   scratch true

   container "${params.gatkContainer}"
   containerOptions "-B ${params.vqsrFolder}:/vqsr"

   input:   
   val flag from vrIndel
   val(vcfFile) from regionsIndels
   
   output:
   val "${vcfFile}" into applyVQSRVCFs

   """
   cp /vqsr/${vcfFile}.excesshet.jvc.vcf.gz* $TMPDIR/
   cp /vqsr/BQC19_indels.recal* $TMPDIR/
   cp /vqsr/BQC19_indels.tranches $TMPDIR/
   cp /vqsr/BQC19_snps.recal* $TMPDIR/
   cp /vqsr/BQC19_snps.tranches $TMPDIR/
   gatk --java-options "-Xmx8g -Xms4g" ApplyVQSR -V $TMPDIR/${vcfFile}.excesshet.jvc.vcf.gz --recal-file $TMPDIR/BQC19_indels.recal --tranches-file $TMPDIR/BQC19_indels.tranches --truth-sensitivity-filter-level 99.7 --create-output-variant-index true -mode INDEL -O $TMPDIR/${vcfFile}.indel.recalibrated.vcf.gz
   gatk --java-options "-Xmx8g -Xms4g" ApplyVQSR -V file:////$TMPDIR/${vcfFile}.indel.recalibrated.vcf.gz --recal-file file:////$TMPDIR/BQC19_snps.recal --tranches-file file:////$TMPDIR/BQC19_snps.tranches --truth-sensitivity-filter-level 99.7 --create-output-variant-index true -mode SNP -O file:////$TMPDIR/${vcfFile}.snp.recalibrated.vcf.gz
   cp $TMPDIR/${vcfFile}.snp.recalibrated.vcf.gz* /vqsr
   """
}*/

/*process ApplyVQSRSNP{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "32 GB" 
   time "4h"
   scratch true

   container "${params.gatkContainer}"
   containerOptions "-B ${params.vqsrFolder}:/vqsr"

   input:   
   val flag from vrSNP
   val(vcfFile) from applyVQSRIndelVCFs
   
   output:
   val "${vcfFile}" into applyVQSRSNPVCFs

   """
   cp /vqsr/${vcfFile}.indel.recalibrated.vcf.gz* $TMPDIR/
   cp /vqsr/BQC19_snps.recal* $TMPDIR/
   cp /vqsr/BQC19_snps.tranches $TMPDIR/
   gatk --java-options "-Xmx8g -Xms4g" ApplyVQSR -V file:////$TMPDIR/${vcfFile}.indel.recalibrated.vcf.gz --recal-file file:////$TMPDIR/BQC19_snps.recal --tranches-file file:////$TMPDIR/BQC19_snps.tranches --truth-sensitivity-filter-level 99.7 --create-output-variant-index true -mode SNP -O file:////$TMPDIR/${vcfFile}.snp.recalibrated.vcf.gz
   cp $TMPDIR/${vcfFile}.snp.recalibrated.vcf.gz* /vqsr
   """
}*/

/*process MergeVQSRVCFs{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "32 GB" 
   time "6h"

   input:   
   val(vcf) from applyVQSRSNPVCFs.collect()

   output:
   val true into result

   """
   bcftools concat -f ${params.vqsrFolder}/vqsrList.txt -a -D -Oz -o ${params.vqsrFolder}/BQC19.vqsr.vcf.gz
   tabix -f -p vcf ${params.vqsrFolder}/BQC19.vqsr.vcf.gz
   """
}*/

