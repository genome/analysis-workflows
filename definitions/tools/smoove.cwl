#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Run Smoove v0.1.6"

baseCommand: "/usr/bin/smoove"
arguments: ["call", "--processes", $(runtime.cores), "-F"]

requirements:
    - class: DockerRequirement
      dockerPull: "quay.io/biocontainers/smoove:0.2.6--0"
    - class: ResourceRequirement
      ramMin: 20000
      coresMin: 4
      tmpdirMin: 10000

inputs:
    bams:
        type: File[]
        inputBinding:
            prefix: --genotype
            position: 1
        doc: "Array of bams to run through lumpy, can be single or small cohort" 
    cohort_name:
        type: string?
        inputBinding:
            prefix: --name
            position: 2
        default: "SV"
        doc: "Used for naming the output file"
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: --fasta
            position: 3
        doc: "Reference used to generate the bams"
    exclude_regions:
        type: File?
        inputBinding:
            prefix: --exclude
            position: 4
        doc: "Bed file used to exclude regions"

outputs:
    output_vcf:
        type: File
        outputBinding:
            glob: "$(inputs.cohort_name)-smoove.genotyped.vcf.gz"
        secondaryFiles: [.csi]
