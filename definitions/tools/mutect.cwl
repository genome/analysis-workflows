#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Mutect2 (GATK 4)"

baseCommand: ["/bin/bash", "Mutect2.sh"]

requirements:
    - class: ResourceRequirement
      ramMin: 32000
      tmpdirMin: 100000
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.2.3.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'Mutect2.sh'
        entry: |
            #!/bin/bash
            set -o pipefail
            set -o errexit

            export tumor_bam="$3"
            export intervals="$4"
            export normal_bam="$5"

            TUMOR=`samtools view -H $tumor_bam | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1`
            if [ "$#" == 4 ]; then
              # tumor only
              /gatk/gatk Mutect2 --java-options "-Xmx20g" -O $1 -R $2 -I "$tumor_bam" -tumor "$TUMOR" -L "$intervals" #Running Mutect2.
            elif [ "$#" == 5 ]; then
              # tumor and normal
              NORMAL=`samtools view -H $normal_bam | perl -nE 'say $1 if /^\@RG.+\tSM:([ -~]+)/' | head -n 1`
              /gatk/gatk Mutect2 --java-options "-Xmx20g" -O $1 -R $2 -I "$tumor_bam" -tumor "$TUMOR" -I "$normal_bam" -normal "$NORMAL" -L "$intervals" #Running Mutect2.
            else
              echo "error in inputs"
              exit 1
            fi

            /gatk/gatk FilterMutectCalls -R $2 -V mutect.vcf.gz -O mutect.filtered.vcf.gz #Running FilterMutectCalls on the output vcf.

arguments:
    - position: 1
      valueFrom: mutect.vcf.gz

inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
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
            position: 5
        secondaryFiles: [.bai]
    interval_list:
        type: File
        inputBinding:
            position: 4

outputs:
    vcf:
        type: File
        outputBinding:
            glob: "mutect.filtered.vcf.gz"
        secondaryFiles: [.tbi]
