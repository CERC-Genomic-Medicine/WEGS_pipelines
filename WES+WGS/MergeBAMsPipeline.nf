bams = Channel.from(file(params.bams_list_path).readLines()).map { line -> fields = line.split(); [ fields[0], file(fields[1]), file(fields[2])] }

process CheckSampleNameMismatch {
   cache "lenient"
   cpus 1
   memory "2 GB"
   time "1h"
   
   input:
   tuple val(sample), file('bam1.bam'),file('bam2.bam') from bams
   
   output:
   tuple val(sample), file('bam1.bam'),file('bam2.bam') into bams_input

   """
   arr=(\$(samtools view -H bam1.bam | grep '^@RG' | sed "s/.*SM:\\([^\t]*\\).*/\1/g"))
   arr+=(\$(samtools view -H bam2.bam | grep '^@RG' | sed "s/.*SM:\\([^\t]*\\).*/\1/g"))
   arrU=(\$(printf "%s\n" "\${arr[@]}" | sort -u))
   if [[ \${#arrU[@]} -gt 1 ]]; then echo "Sample names SM are not unique, please update the bam headers and retry";exit 1;fi
   """

}
process GenerateStatsBeforeMerge {
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "3h"
   errorStrategy "retry"
   maxRetries 3

   input:
   tuple val(sample), file(bam1), file(bam2) from bams_input

   output:
   file "${sample}*bam.*" into bams_stats
   tuple val(sample), file(bam1), file(bam2) into bams_to_merge


   publishDir "${params.result_folder}", pattern: "${sample}*bam.*"

   """
   samtools depth -H -o ${sample}-wes.dedup.bam.depth ${bam1}
   samtools stats ${bam1} > ${sample}-wes.dedup.bam.stats
   samtools coverage -m -o ${sample}-wes.dedup.bam.cov ${bam1}
   samtools depth -H -o ${sample}-wgs.dedup.bam.depth ${bam2}
   samtools stats ${bam2} > ${sample}-wgs.dedup.bam.stats
   samtools coverage -m -o ${sample}-wgs.dedup.bam.cov ${bam2}
   """
}

process MergeAndMarkDuplicates {
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "6h"
   errorStrategy "retry"
   maxRetries 3

   input:
   tuple val(sample), file(bam1), file(bam2) from bams_to_merge

   output:
   tuple file("${sample}.dedup.bam"),file("${sample}.dedup.bai"),file("${sample}.marked_dup_metrics.txt") into bams_merged_dedup
   tuple val(sample), file("${sample}.dedup.bam") into bam_for_stats

   publishDir "${params.result_folder}", pattern: "${sample}.*"

   """
   java -jar -Xmx8G -XX:ParallelGCThreads=1 $EBROOTPICARD/picard.jar MarkDuplicates CREATE_INDEX=true I=${bam1} I=${bam2} O=${sample}.dedup.bam M=${sample}.marked_dup_metrics.txt SORTING_COLLECTION_SIZE_RATIO=0.05 MAX_RECORDS_IN_RAM=100000 COMPRESSION_LEVEL=3 VALIDATION_STRINGENCY=STRICT
   """
}

process GenerateStatsAfterMerge {
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "3h"
   errorStrategy "retry"
   maxRetries 3

   input:
   tuple val(sample), file(mergedBam) from bam_for_stats

   output:
   tuple file("${sample}.dedup.bam.stats"),file("${sample}.dedup.bam.cov"),file("${sample}.dedup.bam.depth") into bams_merged_stats
   tuple val(sample),file(mergedBam) into bam_to_cram

   publishDir "${params.result_folder}", pattern: "${sample}.*"

   """
   samtools depth -H -o ${sample}.dedup.bam.depth ${mergedBam}
   samtools stats ${mergedBam} > ${sample}.dedup.bam.stats
   samtools coverage -m -o ${sample}.dedup.bam.cov ${mergedBam}
   """
}

process GenerateCramFromBam {
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "3h"
   errorStrategy "retry"
   maxRetries 3

   input:
   tuple val(sample),file(mergedBam) from bam_to_cram

   output:
   tuple val(sample),file("${sample}.dedup.cram") into cram

   publishDir "${params.result_folder}", pattern: "${sample}.dedup.cram"

   """
   samtools view -C -T ${params.ref_genome} ${mergedBam} > ${sample}.dedup.cram
   """
}
