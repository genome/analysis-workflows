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

            /gatk/gatk Mutect2 -O $1 -R $2 -I $3 -tumor $5 -I $4 -normal $6 -L $7
            gunzip mutect.vcf.gz
            grep "^##" mutect.vcf >mutect.final.vcf
            grep "^#CHROM" mutect.vcf | sed "s/\t$5\t/\tTUMOR\t/g; s/\t$6/\tNORMAL/g" >> mutect.final.vcf
            grep -v "^#" mutect.vcf >> mutect.final.vcf
            bgzip mutect.final.vcf
            tabix mutect.final.vcf.gz
            mv mutect.vcf.gz.stats mutect.final.vcf.stats
            /gatk/gatk FilterMutectCalls -R $2 -V mutect.final.vcf.gz -O mutect.final.filtered.vcf.gz

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
            position: 7
    tumor_sample_name:
        type: string
        inputBinding:
            position: 5
    normal_sample_name:
        type: string
        inputBinding:
            position: 6

outputs:
    vcf:
        type: File
        outputBinding:
            glob: "mutect.final.filtered.vcf.gz"
        secondaryFiles: [.tbi]
