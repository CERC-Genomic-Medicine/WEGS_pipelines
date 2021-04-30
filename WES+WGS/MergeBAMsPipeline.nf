bams = Channel.from(file(params.bams_list_path).readLines()).map { line -> fields = line.split(); [ fields[0], file(fields[1]), file(fields[2])] }

process CheckSampleNameMismatch {
   cache "lenient"
   cpus 1
   memory "4 GB"
   time "1h"
   
   input:
   tuple val(sample), file('wes_bam.bam'),file('wgs_bam.bam') from bams
   
   output:
   tuple val(sample), file('bam1.bam'), file('bam2.bam') into bams_to_merge
   tuple val(sample), file('bam1.bam') into wes_bam
   tuple val(sample), file('bam1.bam') into wes_before_merge_stats
   tuple val(sample), file('bam2.bam') into (wgs_bam_autosomal, wgs_before_merge_stats)

   """
   samtools view -H wes_bam.bam | sed "s/SM:[^\t]*/SM:${sample}/g" | samtools reheader - wes_bam.bam > bam1.bam
   samtools view -H wgs_bam.bam | sed "s/SM:[^\t]*/SM:${sample}/g" | samtools reheader - wgs_bam.bam > bam2.bam
   """
}

wes_bam.combine(Channel.from(file(params.beds_list_path).readLines()).map { line -> fields = line.split(); [ fields[0], file(fields[1])] }).into { wes_bam_target; wes_bam_flanking }


process MergeAndMarkDuplicates {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "6h"
   errorStrategy "finish"
   

   input:
   tuple val(sample), file(bam1), file(bam2) from bams_to_merge

   output:
   tuple file("${sample}.dedup.bam"),file("${sample}.dedup.bai"),file("${sample}.marked_dup_metrics.txt") into bams_merged_dedup
   tuple val(sample), file("${sample}.dedup.bam") into merged_bam
   tuple val(sample), file("${sample}.dedup.bam") into after_merge_stats
   tuple val(sample), file("${sample}.dedup.bam") into bam_to_cram


   publishDir "${params.result_folder}", pattern: "${sample}.marked_dup_metrics.txt"

   """
   java -jar -Xmx8G -XX:ParallelGCThreads=1 $EBROOTPICARD/picard.jar MarkDuplicates CREATE_INDEX=true I=${bam1} I=${bam2} O=${sample}.dedup.bam M=${sample}.marked_dup_metrics.txt SORTING_COLLECTION_SIZE_RATIO=0.05 MAX_RECORDS_IN_RAM=100000 COMPRESSION_LEVEL=3 VALIDATION_STRINGENCY=STRICT
   """
}

merged_bam.combine(Channel.from(file(params.beds_list_path).readLines()).map { line -> fields = line.split(); [ fields[0], file(fields[1])] }).into { merged_bam_target; merged_bam_flanking }

process GenerateCramFromBam {
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "4h"
   errorStrategy "finish"

   input:
   tuple val(sample),file(mergedBam) from bam_to_cram

   output:
   tuple val(sample),file("${sample}.dedup.cram") into cram_for_indexing

   publishDir "${params.result_folder}", pattern: "${sample}.dedup.cram"

   """
   samtools view -C -T ${params.ref_genome} ${mergedBam} > ${sample}.dedup.cram
   """
}

process IndexCrams {
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "2h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(cram) from cram_for_indexing

   output:
   file("${sample}.dedup.cram.crai") into cram_indexes

   publishDir "${params.result_folder}", pattern: "${sample}.dedup.cram.crai"

   """
   samtools index ${cram} ${sample}.dedup.cram.crai
   """
}

process GenerateStatsBeforeMergeWes {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(bam1) from wes_before_merge_stats

   output:
   file "${sample}-wes.stats" into wes_bef_merge_stats

   publishDir "${params.result_folder}", pattern: "${sample}-wes.stats"

   """
   samtools stats ${bam1} > ${sample}-wes.stats
   """
}

process GenerateStatsBeforeMergeWgs {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(bam2) from wgs_before_merge_stats

   output:
   file "${sample}-wgs.stats" into wgs_bef_merge_stats

   publishDir "${params.result_folder}", pattern: "${sample}-wgs.stats"

   """
   samtools stats ${bam2} > ${sample}-wgs.stats
   """
}

process GenerateStatsAfterMerge {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(mergedBam) from after_merge_stats

   output:
   file "${sample}-merged.stats" into aft_merge_stats

   publishDir "${params.result_folder}", pattern: "${sample}-merged.stats"

   """
   samtools stats ${mergedBam} > ${sample}-merged.stats
   """
}

process GenerateCustomStatsFlankingWes {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(bam1), val(bedLabel), file(bed) from wes_bam_flanking

   output:
   file "${sample}-wes.${bedLabel}.flanking.depth" into wes_flanking_depth

   publishDir "${params.result_folder}", pattern: "${sample}-wes.${bedLabel}.flanking.depth"

   """
   bedtools merge -i ${bed} > ${bedLabel}.adjmerge.bed
   bedtools flank -i ${bedLabel}.adjmerge.bed -g ${params.bed_genome_path} -b 10 > ${bedLabel}.adjmerge.flanking.bed
   samtools depth -a -b ${bedLabel}.adjmerge.flanking.bed -q 20 -Q 20 -s ${bam1} | python ${params.project_dir}/customStats.py -o ${sample}-wes.${bedLabel}.flanking.depth
   """
}

process GenerateCustomStatsTargetWes {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(bam1), val(bedLabel), file(bed) from wes_bam_target

   output:
   file "${sample}-wes.${bedLabel}.depth" into wes_target_depth

   publishDir "${params.result_folder}", pattern: "${sample}-wes.${bedLabel}.depth"

   """
   bedtools merge -i ${bed} > ${bedLabel}.adjmerge.bed
   samtools depth -a -b ${bedLabel}.adjmerge.bed -q 20 -Q 20 -s ${bam1} | python ${params.project_dir}/customStats.py -o ${sample}-wes.${bedLabel}.depth
   """
}
process GenerateCustomStatsFlankingMerged {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(mergedBam), val(bedLabel), file(bed) from merged_bam_flanking

   output:
   file "${sample}-merged.${bedLabel}.flanking.depth" into merged_flanking_depth

   publishDir "${params.result_folder}", pattern: "${sample}-merged.${bedLabel}.flanking.depth"

   """
   bedtools merge -i ${bed} > ${bedLabel}.adjmerge.bed
   bedtools flank -i ${bedLabel}.adjmerge.bed -g ${params.bed_genome_path} -b 10 > ${bedLabel}.adjmerge.flanking.bed
   samtools depth -a -b ${bedLabel}.adjmerge.flanking.bed -q 20 -Q 20 -s ${mergedBam} | python ${params.project_dir}/customStats.py -o ${sample}-merged.${bedLabel}.flanking.depth
   """
}

process GenerateCustomStatsTargetMerged {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(mergedBam), val(bedLabel), file(bed) from merged_bam_target

   output:
   file "${sample}-merged.${bedLabel}.depth" into merged_target_depth

   publishDir "${params.result_folder}", pattern: "${sample}-merged.${bedLabel}.depth"


   """
   bedtools merge -i ${bed} > ${bedLabel}.adjmerge.bed
   samtools depth -a -b ${bedLabel}.adjmerge.bed -q 20 -Q 20 -s ${mergedBam} | python ${params.project_dir}/customStats.py -o ${sample}-merged.${bedLabel}.depth
   """
}

process GenerateCustomStatsWgs {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "3h"
   errorStrategy "finish"

   input:
   tuple val(sample), file(bam2) from wgs_bam_autosomal

   output:
   file "${sample}-wgs.autosomal.depth" into wgs_autosomal_depth

   publishDir "${params.result_folder}", pattern: "${sample}-wgs.autosomal.depth"

   """
   samtools depth -a -q 20 -Q 20 -s ${bam2} | python ${params.project_dir}/customStats.py -a -o ${sample}-wgs.autosomal.depth
   """
}

