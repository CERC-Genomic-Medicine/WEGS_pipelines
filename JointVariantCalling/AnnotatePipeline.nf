regions = Channel.fromList( ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX'] )
//regions1 = Channel.fromList( ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX'] )
//regions2 = Channel.fromList( ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX'] )
//regions = Channel.fromList( ['chr2'] )
//regions = Channel.fromList( ['chr21_1'] )
/*process SplitVCF{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "5 GB" 
   time "10h"

   input:   
   val(chr) from regions

   """
   bcftools view -r ${chr} ${params.multisampleVCF}/BQC19.vqsr.vcf.gz -Oz -o ${params.chrVCF}/BQC19.${chr}.vqsr.vcf.gz
   tabix -f -p vcf ${params.chrVCF}/BQC19.${chr}.vqsr.vcf.gz
   """   

}*/

/*process LeftAlignAndTrim {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "64 GB" 
   time "6h"
   scratch true

   container "${params.gatkContainer}" 
   containerOptions "-B ${params.chrVCF}:/vcfs -B ${params.referenceDir}:/ref -B ${params.resultFolder}:/result"

   input:   
   val(chr) from regions
   
   output:
   val "${chr}" into leftAligned

   """
   cp /vcfs/${params.label}.${chr}.vqsr.vcf.gz* $TMPDIR/
   gatk LeftAlignAndTrimVariants -OVI true -R /ref/${params.referenceGenome}.fa -V $TMPDIR/${params.label}.${chr}.vqsr.vcf.gz -O file:////$TMPDIR/${params.label}.${chr}.leftaligned.vqsr.vcf.gz
   cp $TMPDIR/${params.label}.${chr}.leftaligned.vqsr.vcf.gz* /result
   """
}*/
//zcat $TMPDIR/${params.label}.${chr}.leftaligned.vqsr.vcf.gz | vcf-annotate -a ${params.gatkBundle}/All.renamed_chroms.vcf.gz -c CHROM,POS,ID,REF,ALT,QUAL,FILTER,- --fill-type | bgzip > $TMPDIR/${params.label}.${chr}.annotated.leftaligned.vqsr.vcf.gz
/*process Annotate {

   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "6h"
   scratch true

   input:
   val(chr) from regions

   output:
   val "${chr}" into annotated

   """
   bcftools annotate -a ${params.gatkBundle}/All.renamed_chroms.vcf.gz -c CHROM,POS,ID -Oz -o ${params.resultFolder}/${params.label}.${chr}.annotated.leftaligned.vqsr.vcf.gz ${params.resultFolder}/${params.label}.${chr}.leftaligned.vqsr.vcf.gz
   tabix -f -p vcf ${params.resultFolder}/${params.label}.${chr}.annotated.leftaligned.vqsr.vcf.gz
   """ 

}*/

process SitesOnly {

   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "6h"

   input:
   val(chr) from regions

   output:
   val "${chr}" into sitesonly

   """
   bcftools view -G ${params.resultFolder}/${params.label}.${chr}.indsRemoved.missingFiltered.annotated.leftaligned.vqsr.vcf.gz -Oz -o ${params.resultFolder}/${params.label}.${chr}.sites.indsRemoved.vcf.gz
   tabix -f -p vcf ${params.resultFolder}/${params.label}.${chr}.sites.indsRemoved.vcf.gz
   """ 

}


process MergeSitesOnlyVCFs{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB" 
   time "12h"

   input:   
   val(vcf) from sitesonly.collect()

   output:
   val true into result

   """
   bcftools concat -f ${params.projectFolder}/annotatedVcfList.txt -a -D -Oz -o ${params.chrVCF}/BQC19.sites.indsRemoved.vcf.gz
   tabix -f -p vcf ${params.chrVCF}/BQC19.sites.indsRemoved.vcf.gz
   """
}


/*process VcfStats{

   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"
   scratch true

   input:
   val(chr) from regions

   """
   bcftools filter -e'count(INFO/AC)>1' ${params.resultFolder}/${params.label}.${chr}.annotated.leftaligned.vqsr.vcf.gz | bcftools view -i'FILTER="PASS"' -Oz -o ${params.resultFolder}/${params.label}.${chr}.pass.annotated.leftaligned.vqsr.vcf.gz
   bcftools stats -s- ${params.resultFolder}/${params.label}.${chr}.pass.annotated.leftaligned.vqsr.vcf.gz > ${params.resultFolder}/${params.label}.${chr}.stats
   """ 

} 

process VcfSingletonStats{

   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"
   scratch true

   input:
   val(chr) from regions1

   """
   bcftools filter -e'count(INFO/AC)>1' ${params.resultFolder}/${params.label}.${chr}.annotated.leftaligned.vqsr.vcf.gz | bcftools view -i'FILTER="PASS" && INFO/AC[0]==1' -Oz -o ${params.resultFolder}/${params.label}.${chr}.singleton.annotated.leftaligned.vqsr.vcf.gz
   bcftools stats -s- ${params.resultFolder}/${params.label}.${chr}.singleton.annotated.leftaligned.vqsr.vcf.gz > ${params.resultFolder}/${params.label}.${chr}.singleton.stats
   """ 

} 


process VcfNovelStats{

   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"
   scratch true

   input:
   val(chr) from regions2

   """
   bcftools filter -e'count(INFO/AC)>1' ${params.resultFolder}/${params.label}.${chr}.annotated.leftaligned.vqsr.vcf.gz | bcftools view -i'FILTER="PASS" && ID=="."' -Oz -o ${params.resultFolder}/${params.label}.${chr}.novel.annotated.leftaligned.vqsr.vcf.gz
   bcftools stats -s- ${params.resultFolder}/${params.label}.${chr}.novel.annotated.leftaligned.vqsr.vcf.gz > ${params.resultFolder}/${params.label}.${chr}.novel.stats
   """ 

} */

/*process RemoveFormatFields{

   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "6h"
   scratch true

   input:
   val(chr) from regions

   output:
   val "${chr}" into formatFieldsRemoved

   """
   bcftools annotate -x ^FORMAT/GT,FORMAT/DP ${params.resultFolder}/${params.label}.${chr}.annotated.leftaligned.vqsr.vcf.gz | bcftools view -e'ALT="*"' | bcftools view -i'INFO/AN[0]>1768' -a -s ^BQC12043,BQC12517  -Oz -o ${params.resultFolder}/${params.label}.${chr}.missingFiltered.annotated.leftaligned.vqsr.vcf.gz
   tabix -f -p vcf ${params.resultFolder}/${params.label}.${chr}.missingFiltered.annotated.leftaligned.vqsr.vcf.gz
   """ 

}*/


/*process MissingPerIndv{

   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "6h"
   scratch true

   input:
   val(chr) from regions
   
   output:
   file("${chr}.imiss") into outFiles

   publishDir "${params.resultFolder}"

   """
   vcftools --gzvcf ${params.resultFolder}/${params.label}.${chr}.missingFiltered.annotated.leftaligned.vqsr.vcf.gz --missing-indv --out ${chr}
   """ 

}*/


/*process RemoveIndv{

   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "6h"

   input:
   val(chr) from formatFieldsRemoved
   
   """
   vcftools --gzvcf ${params.resultFolder}/${params.label}.${chr}.missingFiltered.annotated.leftaligned.vqsr.vcf.gz --recode --recode-INFO-all --remove-indv BQC12043 --remove-indv BQC12517 --out ${chr}.indvRemoved
   vcftools --gzvcf ${params.resultFolder}/${params.label}.${chr}.missingFiltered.annotated.leftaligned.vqsr.vcf.gz --singletons --indv BQC12043 --indv BQC12517 --out ${chr}
   cat ${chr}.singletons | grep S | awk -F'\t' '{print \$1,\$2}'|awk 'NR!=1 {print}' > ${chr}.exclude.pos
   vcftools --vcf ${chr}.indvRemoved.recode.vcf  --recode --recode-INFO-all --exclude-positions ${chr}.exclude.pos --out ${chr}.excludePos
   bgzip -c ${chr}.excludePos.recode.vcf > ${params.resultFolder}/${params.label}.${chr}.biallelic.final.vcf.gz
   tabix -p vcf ${params.resultFolder}/${params.label}.${chr}.biallelic.final.vcf.gz
   """ 

} */



/*process CollectVariantCallingMetrics{

   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "6h"

   input:
   val(chr) from regions

   """
   java -jar -Xmx12G $EBROOTPICARD/picard.jar CollectVariantCallingMetrics INPUT=${params.resultFolder}/${params.label}.${chr}.missingFiltered.annotated.leftaligned.vqsr.vcf.gz OUTPUT=${params.chrVCF}/BQC19.${chr} DBSNP=${params.gatkBundle}/All.renamed_chroms.vcf.gz
   """ 

}*/

/*process RemoveMinDP{

   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "12h"

   input:
   val(chr) from regions

   """
   vcftools --gzvcf ${params.resultFolder}/${params.label}.${chr}.missingFiltered.annotated.leftaligned.vqsr.vcf.gz --recode --recode-INFO-all --minDP 30 --out ${chr}.minDP
   vcffixup ${chr}.minDP.recode.vcf > ${chr}.minDP.updated.vcf
   bgzip -c ${chr}.minDP.updated.vcf > ${params.resultFolder}/${params.label}.${chr}.minDP30.vcf.gz
   tabix -p vcf ${params.resultFolder}/${params.label}.${chr}.minDP30.vcf.gz
   """ 

}*/

/*process Annotate {

   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "6h"

   container "${params.gatkContainer}"
   containerOptions "-B ${params.chrVCF}:/out -B ${params.gatkBundle}:/bundle"

   """
   gatk --java-options "-Xmx8g -Xms4g" VariantAnnotator -V /out/BQC19.sites.final.vcf.gz -O /out/BQC19.sites.gatk.vcf.gz --dbsnp /bundle/All.renamed_chroms.vcf.gz
   tabix -f -p vcf /out/BQC19.sites.gatk.vcf.gz
   """ 

}*/


/*process RemoveInd{

   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "6h"
   scratch true

   input:
   val(chr) from regions

   """
   bcftools view -a -s ^BQC10097,BQC10134,BQC10279,BQC10353,BQC10437,BQC10453,BQC10494,BQC10555,BQC10610,BQC10644,BQC10697,BQC10740,BQC10773,BQC10775,BQC10778,BQC10838,BQC10861,BQC10872,BQC10886,BQC10948,BQC11103,BQC11171,BQC11204,BQC11209,BQC11297,BQC11496,BQC11511,BQC11512,BQC11579,BQC11656,BQC11763,BQC11788,BQC11831,BQC11835,BQC11855,BQC11903,BQC11976,BQC12002,BQC12087,BQC12092,BQC12281,BQC12304,BQC12396,BQC12495,BQC12584,BQC12682,BQC12805,BQC12808 ${params.resultFolder}/${params.label}.${chr}.missingFiltered.annotated.leftaligned.vqsr.vcf.gz -Oz -o ${params.resultFolder}/${params.label}.${chr}.indsRemoved.missingFiltered.annotated.leftaligned.vqsr.vcf.gz
   tabix -f -p vcf ${params.resultFolder}/${params.label}.${chr}.indsRemoved.missingFiltered.annotated.leftaligned.vqsr.vcf.gz
   """ 

}*/

//cp ${params.gatkBundle}/All.renamed_chroms.vcf.gz* $TMPDIR/
//cp ${params.resultFolder}/${params.label}.${chr}.leftaligned.vqsr.vcf.gz* $TMPDIR/