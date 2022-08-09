#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2022
*/

inputFiles = Channel.empty()

if( params.inputFileType == "bam" ) {
inputFiles = Channel.fromPath(params.inputFiles).map { file -> [file.getSimpleName(), file.getBaseName(), file, file + ".bai"] }
}
else if( params.inputFileType == "cram" ) {
inputFiles = Channel.fromPath(params.inputFiles).map { file -> [file.getSimpleName(), file.getBaseName(), file, file + ".crai"] }
}

process CramToBam {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "24h"
   
   input:   
   tuple val(input_label), val(filename), file(inputFile), file(index) from inputFiles

   output:
   tuple val(input_label), val(filename) into (for_Delly, for_manta, for_BD, for_BS, for_CN, for_Lumpy, melt_alu, melt_hervk, melt_line1, melt_sva)
   
   script:
   if (params.doCramToBam)
      """
      module load samtools
      
      mkdir -p ${params.bamsFolder}
      samtools view -b -T ${params.referenceGenome} ${inputFile} > ${params.bamsFolder}/${filename}.bam
      samtools index -b ${params.bamsFolder}/${filename}.bam
      """
   else
      """
      echo ${input_label}
      """

}

process DellySV {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "64 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), val(filename) from for_Delly

   """
   module load bcftools
   module load samtools
   
   mkdir -p ${params.result}/delly
   export OMP_NUM_THREADS=1
   ${params.delly} call -g ${params.referenceGenome} -o ${input_label}.bcf ${params.bamsFolder}/${filename}.bam
   bcftools view ${input_label}.bcf -Oz -o ${input_label}.vcf.gz
   tabix -f -p vcf ${input_label}.vcf.gz
   cp ${input_label}.vcf.gz* ${params.result}/delly
   """
}

process MantaSV {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), val(filename) from for_manta

   """
   module load gcc
   module load manta
   
   mkdir work
   mkdir -p ${params.result}/manta
   ${params.manta}/configManta.py --bam ${params.bamsFolder}/${filename}.bam --referenceFasta ${params.referenceGenome} --callRegions ${params.referenceDir}/GRCh38.include.bed.gz  --runDir work
   ./work/runWorkflow.py
   mv ./work/results/variants/diploidSV.vcf.gz ${params.result}/manta/${input_label}.vcf.gz
   mv ./work/results/variants/diploidSV.vcf.gz.tbi ${params.result}/manta/${input_label}.vcf.gz.tbi
   """

}


process BreakDancerSV {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), val(filename) from for_BD

   """
   module load python
   module load bcftools
   export PYTHONPATH="\${PYTHONPATH}:${params.sve}/:${params.sve}/bin/:${params.sve}/stages/"
   
   mkdir work
   mkdir -p ${params.result}/breakdancer
   ${params.sve}/bin/sve call -r ${params.referenceGenome} -g ${params.genome} -a breakdancer -o work ${params.bamsFolder}/${filename}.bam
   bcftools view ./work/${filename}_S4.vcf -Oz -o ./work/${input_label}.vcf.gz
   tabix -f -p vcf ./work/${input_label}.vcf.gz
   mkdir -p ${params.result}/breakdancer
   cp -r ./work/${input_label}.vcf.gz* ${params.result}/breakdancer
   """

}
process BreakSeqSV {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), val(filename) from for_BS

   """
   module load python
   module load bcftools
   export PYTHONPATH="\${PYTHONPATH}:${params.sve}/:${params.sve}/bin/:${params.sve}/stages/"   
   mkdir work
   mkdir -p ${params.result}/breakseq
   ${params.sve}/bin/sve call -r ${params.referenceGenome} -g ${params.genome} -a breakseq -o work ${params.bamsFolder}/${filename}.bam
   bcftools view ./work/${filename}_S35.vcf -Oz -o ./work/${input_label}.vcf.gz
   tabix -f -p vcf ./work/${input_label}.vcf.gz
   mkdir -p ${params.result}/breakseq
   cp -r ./work/${input_label}.vcf.gz* ${params.result}/breakseq
   """

}
process CNVnatorSV {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "32 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), val(filename) from for_CN

   """
   module load python
   module load bcftools
   export PYTHONPATH="\${PYTHONPATH}:${params.sve}/:${params.sve}/bin/:${params.sve}/stages/"   
   mkdir work
   mkdir -p ${params.result}/cnvnator
   ${params.sve}/bin/sve call -r ${params.referenceGenome} -M 24 -g ${params.genome} -a cnvnator -o work ${params.bamsFolder}/${filename}.bam
   bcftools view ./work/${filename}_S10.vcf -Oz -o ./work/${input_label}.vcf.gz
   tabix -f -p vcf ./work/${input_label}.vcf.gz
   mkdir -p ${params.result}/cnvnator
   cp -r ./work/${input_label}.vcf.gz* ${params.result}/cnvnator
   """

}

process LumpySV {
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "64 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), val(filename) from for_Lumpy

   """
   module load samtools
   module load bcftools
   
   mkdir work
   mkdir -p ${params.result}/lumpy
   samtools view -b -F 1294 ${params.bamsFolder}/${filename}.bam > ./work/${input_label}.discordants.bam
   samtools view -h ${params.bamsFolder}/${filename}.bam | ${params.lumpy}/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > ${input_label}.splitters.bam

   ${params.lumpy}/bin/lumpyexpress -B ${params.bamsFolder}/${filename}.bam -S ./work/${input_label}.discordants.bam -D ./${input_label}.splitters.bam -o ./work/${input_label}.vcf
   cat ./work/${input_label}.vcf | awk '\$1 ~ /^#/ {print \$0;next} {print \$0 | "sort -k1,1 -k2,2n"}' > ./work/${input_label}.sorted.vcf
   bgzip < ./work/${input_label}.sorted.vcf > ./work/${input_label}.vcf.gz 
   tabix -f -p vcf ./work/${input_label}.vcf.gz
   mkdir -p ${params.result}/lumpy
   cp -r ./work/${input_label}.vcf.gz* ${params.result}/lumpy
   """

}

process MeltALU{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "96h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), val(filename) from melt_alu

   """
   module load bowtie2
   module load java
   module load samtools
   
   mkdir work
   mkdir -p ./work/ALU
   mkdir -p ${params.result}/melt/ALU   
   cp ${params.bamsFolder}/${filename}.bam ./work
   cp ${params.bamsFolder}/${filename}*.bai ./work
   java -Xmx8G -jar ${params.melt}/MELT.jar Preprocess -bamfile ./work/${filename}.bam -h ${params.referenceGenome}
   java -Xmx12G -jar ${params.melt}/MELT.jar IndivAnalysis -bamfile ./work/${filename}.bam -w ./work/ALU -t ${params.mei}/ALU_MELT.zip -h ${params.referenceGenome}
   mkdir -p ${params.result}/melt/ALU
   rm -r ./work/ALU/*tmp
   cp -r ./work/ALU/* ${params.result}/melt/ALU
   """

}

process MeltHERVK{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "96h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), val(filename) from melt_hervk

   """
   module load bowtie2
   module load java
   module load samtools
   
   mkdir work
   mkdir -p ./work/HERVK
   mkdir -p ${params.result}/melt/HERVK
   cp ${params.bamsFolder}/${filename}.bam ./work
   cp ${params.bamsFolder}/${filename}*.bai ./work
   java -Xmx8G -jar ${params.melt}/MELT.jar Preprocess -bamfile ./work/${filename}.bam -h ${params.referenceGenome}
   java -Xmx12G -jar ${params.melt}/MELT.jar IndivAnalysis -bamfile ./work/${filename}.bam -w ./work/HERVK -t ${params.mei}/HERVK_MELT.zip -h ${params.referenceGenome}
   mkdir -p ${params.result}/melt/HERVK
   rm -r ./work/HERVK/*tmp
   cp -r ./work/HERVK/* ${params.result}/melt/HERVK
   """

}

process MeltLINE1{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "96h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), val(filename) from melt_line1

   """
   module load bowtie2
   module load java
   module load samtools
   
   mkdir work
   mkdir -p ./work/LINE1
   mkdir -p ${params.result}/melt/LINE1
   
   cp ${params.bamsFolder}/${filename}.bam ./work
   cp ${params.bamsFolder}/${filename}*.bai ./work
   java -Xmx8G -jar ${params.melt}/MELT.jar Preprocess -bamfile ./work/${filename}.bam -h ${params.referenceGenome}
   java -Xmx12G -jar ${params.melt}/MELT.jar IndivAnalysis -bamfile ./work/${filename}.bam -w ./work/LINE1 -t ${params.mei}/LINE1_MELT.zip -h ${params.referenceGenome}
   mkdir -p ${params.result}/melt/LINE1
   rm -r ./work/LINE1/*tmp
   cp -r ./work/LINE1/* ${params.result}/melt/LINE1
   """

}

process MeltSVA{
   errorStrategy 'retry'
   maxRetries 3
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "96h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), val(filename) from melt_sva

   """
   module load bowtie2
   module load java
   module load samtools
   
   mkdir work
   mkdir -p ./work/SVA
   mkdir -p ${params.result}/melt/SVA
   
   cp ${params.bamsFolder}/${filename}.bam ./work
   cp ${params.bamsFolder}/${filename}*.bai ./work
   java -Xmx8G -jar ${params.melt}/MELT.jar Preprocess -bamfile ./work/${filename}.bam -h ${params.referenceGenome}
   java -Xmx12G -jar ${params.melt}/MELT.jar IndivAnalysis -bamfile ./work/${filename}.bam -w ./work/SVA -t ${params.mei}/SVA_MELT.zip -h ${params.referenceGenome}
   mkdir -p ${params.result}/melt/SVA
   rm -r ./work/SVA/*tmp
   cp -r ./work/SVA/* ${params.result}/melt/SVA
   """

}