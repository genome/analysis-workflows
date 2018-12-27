#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "add bam_readcount info to vcf"

baseCommand: ["vcf-readcount-annotator"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/vcf_annotation_tools-cwl:1.4.6"
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
        type: string
        inputBinding:
            position: 3
    sample_name:
        type: string
        inputBinding:
            prefix: "-s"
outputs:
    annotated_bam_readcount_vcf:
        type: File
        outputBinding:
            glob: "annotated.bam_readcount.vcf.gz"

