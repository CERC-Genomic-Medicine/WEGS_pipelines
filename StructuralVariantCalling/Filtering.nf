#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2022
*/

process HarmonizeMELT {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "6h"

   output:
   val 1 into merge

   """
   module load r
   module load vcftools
   module load bcftools
   module load samtools
   module load picard
   module load bedtools

   mkdir -p ${params.result}/filtering/fromMELT
   Rscript ${params.scripts}/harmonizeInsertions.R ${params.result} samples_merged_INS.raw.vcf.gz

   Rscript ${params.scripts}/filterINS.R samples_merged_INS.raw.vcf.gz samples_merged_INS.raw2.vcf.gz
   zcat samples_merged_INS.raw2.vcf.gz | bcftools sort -Oz -o samples_merged_INS.raw.vcf.gz

   Rscript ${params.scripts}/removeSampleOutliers.R samples_merged_INS.raw.vcf.gz samples_merged_INS.filt1.vcf.gz outlier_samples.INS.blacklist


   cat outlier_samples.INS.blacklist | cut -f1 > outlier_ids
   bcftools query -l samples_merged_INS.raw.vcf.gz | grep -v -x -f outlier_ids > samplesToMerge.list
   zcat samples_merged_INS.raw.vcf.gz | bcftools view --force-samples -S samplesToMerge.list -Oz > samples_merged_INS.filt2.vcf.gz

   bcftools view -c1 samples_merged_INS.filt2.vcf.gz | vcftools --vcf - --max-missing 0.10 --recode --recode-INFO-all --stdout | sed '/^##bcftools/d' | bcftools sort -Oz -o samples_merged_INS.filt3.vcf.gz

   Rscript ${params.scripts}/renameSVIds.R samples_merged_INS.filt3.vcf.gz samples_merged_INS.filt4.vcf.gz

   zcat samples_merged_INS.filt4.vcf.gz | bcftools +fill-tags -Oz -o samples_merged_INS.Final.vcf.gz

   bcftools view -q 0.05:minor samples_merged_INS.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o samples_merged_INS.maf05.Final.vcf.gz
   bcftools view -q 0.01:minor samples_merged_INS.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o samples_merged_INS.maf01.Final.vcf.gz

   tabix -f -p vcf samples_merged_INS.maf05.Final.vcf.gz
   tabix -f -p vcf samples_merged_INS.maf01.Final.vcf.gz

   bcftools view samples_merged_INS.Final.vcf.gz | bcftools +missing2ref | bcftools +fill-tags | bcftools sort | bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%SVLEN\t%AC\t%AF\t%AN\n' > INS.dat
   bedtools annotate -both -i INS.dat -files ${params.scripts}/hg38.SegmentalDups.bed.gz ${params.scripts}/hg38.SimpleRepeats.bed.gz > INS.sd_sr_cov.txt

   cp INS.dat  ${params.result}/filtering/fromMELT
   cp INS.sd_sr_cov.txt ${params.result}/filtering/fromMELT
   cp samples_merged_INS.Final.vcf.gz ${params.result}/filtering/fromMELT
   cp samples_merged_INS.maf* ${params.result}/filtering/fromMELT
   """
}

process MergeSmoove {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "2h"

   input:
   val(flag) from merge

   """
   module load r
   module load vcftools
   module load bcftools
   module load samtools
   module load picard
   module load bedtools

   mkdir -p ${params.result}/filtering/fromSMOOVE


   cp ${params.result}/genotype/Results/DEL.smoove.square.vcf.gz  ${params.result}/filtering/fromSMOOVE/samples_merged_DEL.raw.vcf.gz
   cp ${params.result}/genotype/Results/DUP.smoove.square.vcf.gz  ${params.result}/filtering/fromSMOOVE/samples_merged_DUP.raw.vcf.gz
   cp ${params.result}/genotype/Results/INV.smoove.square.vcf.gz  ${params.result}/filtering/fromSMOOVE/samples_merged_INV.raw.vcf.gz
   cp ${params.result}/mergeSamples/samples_merged_TRABND.vcf ${params.result}/filtering/fromSMOOVE/samples_merged_TRABND.vcf

   bgzip -f -c ${params.result}/mergeSamples/samples_merged_TRABND.vcf > ${params.result}/filtering/fromSMOOVE/samples_merged_TRABND.vcf.gz
   tabix -p vcf ${params.result}/filtering/fromSMOOVE/samples_merged_TRABND.vcf.gz

   Rscript ${params.scripts}/removeSampleOutliers.R ${params.result}/filtering/fromSMOOVE/samples_merged_DEL.raw.vcf.gz samples_merged_DEL.filt1.vcf.gz outlier_samples.DEL.SM.blacklist
   Rscript ${params.scripts}/removeSampleOutliers.R ${params.result}/filtering/fromSMOOVE/samples_merged_DUP.raw.vcf.gz samples_merged_DUP.filt1.vcf.gz outlier_samples.DUP.SM.blacklist
   Rscript ${params.scripts}/removeSampleOutliers.R ${params.result}/filtering/fromSMOOVE/samples_merged_INV.raw.vcf.gz samples_merged_INV.filt1.vcf.gz outlier_samples.INV.SM.blacklist

   cat outlier_samples.DEL.SM.blacklist outlier_samples.DUP.SM.blacklist outlier_samples.INV.SM.blacklist | sort | uniq > outlier_samples.blacklist
   cat outlier_samples.blacklist | cut -f1 > outlier_ids
   bcftools query -l ${params.result}/filtering/fromSMOOVE/samples_merged_DEL.raw.vcf.gz | grep -v -x -f outlier_ids > samplesToMerge.list
   
   bcftools view --force-samples -S samplesToMerge.list ${params.result}/filtering/fromSMOOVE/samples_merged_DEL.raw.vcf.gz -Oz > samples_merged_DEL.filt2.vcf.gz
   bcftools view --force-samples -S samplesToMerge.list ${params.result}/filtering/fromSMOOVE/samples_merged_DUP.raw.vcf.gz -Oz > samples_merged_DUP.filt2.vcf.gz
   bcftools view --force-samples -S samplesToMerge.list ${params.result}/filtering/fromSMOOVE/samples_merged_INV.raw.vcf.gz -Oz > samples_merged_INV.filt2.vcf.gz

   bcftools view -c1 samples_merged_DEL.filt2.vcf.gz | vcftools --vcf - --max-missing 0.10 --recode --recode-INFO-all --stdout | sed '/^##bcftools/d' | bcftools sort -Oz -o samples_merged_DEL.filt3.vcf.gz
   bcftools view -c1 samples_merged_DUP.filt2.vcf.gz | vcftools --vcf - --max-missing 0.10 --recode --recode-INFO-all --stdout | sed '/^##bcftools/d' | bcftools sort -Oz -o samples_merged_DUP.filt3.vcf.gz 
   bcftools view -c1 samples_merged_INV.filt2.vcf.gz | vcftools --vcf - --max-missing 0.10 --recode --recode-INFO-all --stdout | sed '/^##bcftools/d' | bcftools sort -Oz -o samples_merged_INV.filt3.vcf.gz 

   Rscript ${params.scripts}/filterSVByFC.R samples_merged_DEL.filt3.vcf.gz samples_merged_DEL.filt4.vcf.gz
   Rscript ${params.scripts}/filterSVByFC.R samples_merged_DUP.filt3.vcf.gz samples_merged_DUP.filt4.vcf.gz
   Rscript ${params.scripts}/filterSVByFC.R samples_merged_INV.filt3.vcf.gz samples_merged_INV.filt4.vcf.gz

   Rscript ${params.scripts}/renameSVIds.R samples_merged_DEL.filt4.vcf.gz samples_merged_DEL.filt5.vcf.gz
   Rscript ${params.scripts}/renameSVIds.R samples_merged_DUP.filt4.vcf.gz samples_merged_DUP.filt5.vcf.gz
   Rscript ${params.scripts}/renameSVIds.R samples_merged_INV.filt4.vcf.gz samples_merged_INV.filt5.vcf.gz

   zcat samples_merged_DEL.filt5.vcf.gz | bcftools annotate -x "INFO/SUPP_VEC,INFO/SUPP" > samples_merged_DEL.filt6.vcf
   zcat samples_merged_DUP.filt5.vcf.gz | bcftools annotate -x "INFO/SUPP_VEC,INFO/SUPP" > samples_merged_DUP.filt6.vcf
   zcat samples_merged_INV.filt5.vcf.gz | bcftools annotate -x "INFO/SUPP_VEC,INFO/SUPP" > samples_merged_INV.filt6.vcf

   java -jar \$EBROOTPICARD/picard.jar FixVcfHeader I=samples_merged_DEL.filt6.vcf O=samples_merged_DEL.filt7.vcf
   java -jar \$EBROOTPICARD/picard.jar FixVcfHeader I=samples_merged_DUP.filt6.vcf O=samples_merged_DUP.filt7.vcf
   java -jar \$EBROOTPICARD/picard.jar FixVcfHeader I=samples_merged_INV.filt6.vcf O=samples_merged_INV.filt7.vcf

   bgzip -f samples_merged_DEL.filt7.vcf
   bgzip -f samples_merged_DUP.filt7.vcf
   bgzip -f samples_merged_INV.filt7.vcf

   bcftools +fill-tags -Oz -o samples_merged_DEL.Final.vcf.gz samples_merged_DEL.filt7.vcf.gz
   bcftools +fill-tags -Oz -o samples_merged_DUP.Final.vcf.gz samples_merged_DUP.filt7.vcf.gz
   bcftools +fill-tags -Oz -o samples_merged_INV.Final.vcf.gz samples_merged_INV.filt7.vcf.gz

   bcftools view -q 0.05:minor samples_merged_DEL.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o samples_merged_DEL.rate90.maf05.Final.vcf.gz
   bcftools view -q 0.05:minor samples_merged_DUP.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o samples_merged_DUP.rate90.maf05.Final.vcf.gz
   bcftools view -q 0.05:minor samples_merged_INV.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o samples_merged_INV.rate90.maf05.Final.vcf.gz

   bcftools view -q 0.01:minor samples_merged_DEL.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o samples_merged_DEL.rate90.maf01.Final.vcf.gz
   bcftools view -q 0.01:minor samples_merged_DUP.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o samples_merged_DUP.rate90.maf01.Final.vcf.gz
   bcftools view -q 0.01:minor samples_merged_INV.Final.vcf.gz | bcftools +fill-tags | bcftools sort -Oz -o samples_merged_INV.rate90.maf01.Final.vcf.gz

   tabix -f -p vcf samples_merged_DEL.rate90.maf05.Final.vcf.gz
   tabix -f -p vcf samples_merged_DUP.rate90.maf05.Final.vcf.gz
   tabix -f -p vcf samples_merged_INV.rate90.maf05.Final.vcf.gz

   tabix -f -p vcf samples_merged_DEL.rate90.maf01.Final.vcf.gz
   tabix -f -p vcf samples_merged_DUP.rate90.maf01.Final.vcf.gz
   tabix -f -p vcf samples_merged_INV.rate90.maf01.Final.vcf.gz

   bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%SVLEN\t%AC\t%AF\t%AN\n' samples_merged_DEL.Final.vcf.gz > DEL.dat
   bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%SVLEN\t%AC\t%AF\t%AN\n' samples_merged_DUP.Final.vcf.gz > DUP.dat
   bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%SVTYPE\t%SVLEN\t%AC\t%AF\t%AN\n' samples_merged_INV.Final.vcf.gz > INV.dat

   bedtools annotate -both -i DEL.dat -files ${params.scripts}/hg38.SegmentalDups.bed.gz ${params.scripts}/hg38.SimpleRepeats.bed.gz > DEL.sd_sr_cov.txt
   bedtools annotate -both -i DUP.dat -files ${params.scripts}/hg38.SegmentalDups.bed.gz ${params.scripts}/hg38.SimpleRepeats.bed.gz > DUP.sd_sr_cov.txt
   bedtools annotate -both -i INV.dat -files ${params.scripts}/hg38.SegmentalDups.bed.gz ${params.scripts}/hg38.SimpleRepeats.bed.gz > INV.sd_sr_cov.txt

   java -jar \$EBROOTPICARD/picard.jar MergeVcfs I=samples_merged_DEL.Final.vcf.gz I=samples_merged_DUP.Final.vcf.gz I=samples_merged_INV.Final.vcf.gz I=${params.result}/filtering/fromMELT/samples_merged_INS.Final.vcf.gz I=${params.result}/filtering/fromSMOOVE/samples_merged_TRABND.vcf.gz O=samples_merged_ALL.Final.vcf.gz
   tabix -p vcf samples_merged_ALL.Final.vcf.gz

   cp *.dat  ${params.result}/filtering/fromSMOOVE
   cp *.sd_sr_cov.txt ${params.result}/filtering/fromSMOOVE
   cp samples_merged_*.Final.vcf.gz* ${params.result}/filtering/fromSMOOVE
   cp samples_merged_*.maf* ${params.result}/filtering/fromSMOOVE
   """
}

