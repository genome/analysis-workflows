#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "add VEP annotation to report"

baseCommand: ["vep-annotation-reporter"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/vcf_annotation_tools-cwl:3.0.0"
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/variants.annotated.tsv }]
inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
    vep_fields:
        type:
            type: array
            items: string
        default: ['Consequence','SYMBOL','Feature','HGVSc','HGVSp']
        inputBinding:
            position: 2
    tsv:
        type: File?
        inputBinding:
            prefix: "-t"
outputs:
    annotated_variants_tsv:
        type: File
        outputBinding:
            glob: variants.annotated.tsv
