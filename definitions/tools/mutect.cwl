#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Mutect2 (GATK 4)"

baseCommand: ["/bin/bash", "Mutect2.sh"]

requirements:
    - class: ResourceRequirement
      ramMin: 20000
      tmpdirMin: 100000
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.2.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'Mutect2.sh'
        entry: |
            set -o pipefail
            set -o errexit
            
            TUMOR=`(samtools view -H $3 | grep "^@RG" | sed "s/.*SM:\([^\t]*\).*/\1/g" | head -n 1)` #Extracting the sample name from the TUMOR bam.
            NORMAL=`(samtools view -H $4 | grep "^@RG" | sed "s/.*SM:\([^\t]*\).*/\1/g" | head -n 1)` #Extracting the sample name from the NORMAL bam.
            /gatk/gatk Mutect2 -O $1 -R $2 -I $3 -tumor $TUMOR -I $4 -normal $NORMAL -L $5 #Running Mutect2.
            gunzip mutect.vcf.gz #Unzipping the Mutect2 output vcf file to alter the header names.
            grep "^##" mutect.vcf >mutect.final.vcf #Extracting all the lines above #CHROM in the vcf to a new output vcf.
            grep "^#CHROM" mutect.vcf | sed "s/\<$TUMOR\>/TUMOR/g; s/\<$NORMAL\>/NORMAL/g" >> mutect.final.vcf #Concatenating the #CHROM line with substituted sample headers to the output vcf.
            grep -v "^#" mutect.vcf >> mutect.final.vcf #Concatenating the sites to the output vcf.
            bgzip mutect.final.vcf #Bzipping the final output vcf.
            tabix mutect.final.vcf.gz #Indexing the bzipped output vcf.
            mv mutect.vcf.gz.stats mutect.final.vcf.gz.stats #Renaming the .stats file to match the name conversion of the final output vcf.
            /gatk/gatk FilterMutectCalls -R $2 -V mutect.final.vcf.gz -O mutect.final.filtered.vcf.gz #Running FilterMutectCalls on the output vcf.

arguments:
    - position: 1
      valueFrom: mutect.vcf.gz

inputs:
    reference:
        type: string
        inputBinding:
            position: 2
    tumor_bam:
        type: File
        inputBinding:
            position: 3
        secondaryFiles: [.bai]
    normal_bam:
        type: File?
        inputBinding:
            position: 4
        secondaryFiles: [.bai]
    interval_list:
        type: File
        inputBinding:
            position: 5

outputs:
    vcf:
        type: File
        outputBinding:
            glob: "mutect.final.filtered.vcf.gz"
        secondaryFiles: [.tbi]
