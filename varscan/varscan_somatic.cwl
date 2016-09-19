#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "varscan v2.4.2 somatic"
baseCommand: "/usr/bin/cwl_helper.sh"
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/varscan-cwl:v2.4.2-samtools1.3.1"
inputs:
    tumor_bam:
        type: File
        inputBinding:
            position: 1
    normal_bam:
        type: File
        inputBinding:
            position: 2
    reference:
        type: File
        inputBinding:
            position: 3
        secondaryFiles: .fai
outputs:
    snvs:
        type: File
        outputBinding:
            glob: "output.snp.vcf"
    indels:
        type: File
        outputBinding:
            glob: "output.indel.vcf"
