#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "varscan v2.4.2 somatic"
baseCommand: "somatic"
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/varscan:v2.4.2"
    - class: InlineJavascriptRequirement
arguments:
    - { valueFrom: $(runtime.outdir + "/output") }
    - "--mpileup"
    - "1"
    - "--output-vcf"
inputs:
    mpileup:
        type: File
        inputBinding:
            position: -1 #before arguments
        streamable: true
outputs:
    snvs:
        type: File
        outputBinding:
            glob: "output.snp.vcf"
    indels:
        type: File
        outputBinding:
            glob: "output.indel.vcf"
