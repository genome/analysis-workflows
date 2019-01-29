#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "add bam_readcount info to vcf"

baseCommand: ["vcf-readcount-annotator"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/vcf_annotation_tools-cwl:3.0.0"
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/annotated.bam_readcount.vcf.gz }]
inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
    bam_readcount_tsv:
        type: File
        inputBinding:
            position: 2
    data_type:
        type:
            - type: enum
              symbols: ["DNA", "RNA"]
        inputBinding:
            position: 3
    variant_type:
        type:
            - "null"
            - type: enum
              symbols: ["snv", "indel", "all"]
        inputBinding:
            prefix: "-t"
    sample_name:
        type: string?
        inputBinding:
            prefix: "-s"
outputs:
    annotated_bam_readcount_vcf:
        type: File
        outputBinding:
            glob: "annotated.bam_readcount.vcf.gz"

