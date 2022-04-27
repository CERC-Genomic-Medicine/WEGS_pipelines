#!/usr/bin/env nextflow

/*
*AUTHOR: Praveen Nadukkalam Ravindran, PhD <praveen.nadukkalamravindran@mcgill.ca>
*VERSION: 1.0
*YEAR: 2022
*/

inputFiles = Channel.fromPath(params.inputFiles).map { file -> [file.getSimpleName(), file, file + ".crai"] }
//inputFiles = Channel.from(file(params.inputFiles).readLines()).map { line -> fields = line.split(); [ fields[0], file(fields[1])] }


/*process DellySV {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "64 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), file(inputFile) from inputFiles

   """
   module load bcftools
   module load samtools
   
   export OMP_NUM_THREADS=1
   samtools index ${inputFile}
   ${params.delly} call -g ${params.referenceGenome} -o ${input_label}.bcf ${inputFile}
   bcftools view ${input_label}.bcf -Oz -o ${input_label}.vcf.gz
   tabix -f -p vcf ${input_label}.vcf.gz
   cp ${input_label}.vcf.gz* ${params.result}/delly
   """
}*/

/*process MantaSV {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "32 GB"
   time "96h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), file(inputFile), file(index) from inputFiles

   """
   module load gcc
   module load manta
   
   mkdir work
   ${params.manta}/configManta.py --bam ${inputFile} --referenceFasta ${params.referenceGenome} --callRegions ${params.referenceDir}/GRCh38.include.bed.gz  --runDir work
   ./work/runWorkflow.py
   mkdir ${params.result}/manta/${input_label}
   cp -r ./work/results/* ${params.result}/manta/${input_label}
   """

}*/

process CramToBam {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "12h"
   //scratch '$SLURM_TMPDIR'
   
   input:   
   tuple val(input_label), file(inputFile), file(index) from inputFiles

   output:
   tuple val(input_label), file("${input_label}.bam"), file("${input_label}.bam.bai") into (for_BD, for_BS, for_CN, for_Lumpy, for_Melt  )


   """
   module load samtools

   samtools view -b -T ${params.referenceGenome} ${inputFile} > ${input_label}.bam
   samtools index -b ${input_label}.bam
   """

}


process BreakDancerSV {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), file(inputFile), file(index) from for_BD

   """
   module load python
   module load bcftools
   export PYTHONPATH="\${PYTHONPATH}:/lustre03/project/rrg-vmooser/praveen/tools/SVE/:/lustre03/project/rrg-vmooser/praveen/tools/SVE/bin/:/lustre03/project/rrg-vmooser/praveen/tools/SVE/stages/"
   
   mkdir work
   ${params.sve}/sve call -r ${params.referenceGenome} -g ${params.genome} -a breakdancer -o work ${inputFile}
   bcftools view ./work/${input_label}_S4.vcf -Oz -o ./work/${input_label}.vcf.gz
   tabix -f -p vcf ./work/${input_label}.vcf.gz
   mkdir -p ${params.result}/breakdancer
   cp -r ./work/${input_label}.vcf.gz* ${params.result}/breakdancer
   """

}
process BreakSeqSV {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), file(inputFile), file(index) from for_BS

   """
   module load python
   module load bcftools
   export PYTHONPATH="\${PYTHONPATH}:/lustre03/project/rrg-vmooser/praveen/tools/SVE/:/lustre03/project/rrg-vmooser/praveen/tools/SVE/bin/:/lustre03/project/rrg-vmooser/praveen/tools/SVE/stages/"
   
   mkdir work
   ${params.sve}/sve call -r ${params.referenceGenome} -g ${params.genome} -a breakseq -o work ${inputFile}
   bcftools view ./work/${input_label}_S35.vcf -Oz -o ./work/${input_label}.vcf.gz
   tabix -f -p vcf ./work/${input_label}.vcf.gz
   mkdir -p ${params.result}/breakseq
   cp -r ./work/${input_label}.vcf.gz* ${params.result}/breakseq
   """

}
process CNVnatorSV {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "32 GB"
   time "72h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), file(inputFile), file(index) from for_CN

   """
   module load python
   module load bcftools
   export PYTHONPATH="\${PYTHONPATH}:/lustre03/project/rrg-vmooser/praveen/tools/SVE/:/lustre03/project/rrg-vmooser/praveen/tools/SVE/bin/:/lustre03/project/rrg-vmooser/praveen/tools/SVE/stages/"
   
   mkdir work
   ${params.sve}/sve call -r ${params.referenceGenome} -M 24 -g ${params.genome} -a cnvnator -o work ${inputFile}
   bcftools view ./work/${input_label}_S10.vcf -Oz -o ./work/${input_label}.vcf.gz
   tabix -f -p vcf ./work/${input_label}.vcf.gz
   mkdir -p ${params.result}/cnvnator
   cp -r ./work/${input_label}.vcf.gz* ${params.result}/cnvnator
   """

}

process LumpySV {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "16 GB"
   time "48h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), file(inputFile), file(index) from for_Lumpy

   """
   module load samtools
   module load bcftools
   
   mkdir work
   samtools view -b -F 1294 ${inputFile} > ./work/${input_label}.discordants.bam
   samtools view -h ${inputFile} | ${params.lumpy}/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > ${input_label}.splitters.bam

   ${params.lumpy}/bin/lumpyexpress -B ${inputFile} -S ./work/${input_label}.discordants.bam -D ./${input_label}.splitters.bam -o ./work/${input_label}.vcf
   mkdir -p ${params.result}/lumpy
   cp -r ./work/${input_label}.vcf ${params.result}/lumpy
   """

}

/*process MeltSV {
   errorStrategy 'ignore'
   cache "lenient"
   cpus 1
   memory "8 GB"
   time "2h"
   scratch '$SLURM_TMPDIR'

   input:   
   tuple val(input_label), file(inputFile), file(index) from for_Melt

   """
   module load samtools
   module load bcftools
   
   mkdir work
   samtools view -b -F 1294 ${inputFile} > ./work/${input_label}.discordants.bam
   samtools view -h ${inputFile} | ${params.lumpy}/scripts/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > ${input_label}.splitters.bam

   ${params.lumpy}/bin/lumpyexpress -B ${inputFile} -S ./work/${input_label}.discordants.bam -D ./${input_label}.splitters.bam -o ./work/${input_label}.vcf
   bcftools view ./work/${input_label}.vcf -Oz -o ./work/${input_label}.vcf.gz
   tabix -f -p vcf ./work/${input_label}.vcf.gz
   mkdir -p ${params.result}/lumpy
   cp -r ./work/${input_label}.vcf.gz* ${params.testResult}/lumpy
   """

}*/

