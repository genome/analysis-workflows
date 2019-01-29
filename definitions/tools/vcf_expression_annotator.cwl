#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "add expression info to vcf"

baseCommand: ["vcf-expression-annotator"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/vcf_annotation_tools-cwl:3.0.0"
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/annotated.expression.vcf.gz }]
inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
    expression_file:
        type: File
        inputBinding:
            position: 2
    expression_tool:
        type: string
        inputBinding:
            position: 3
    data_type:
        type: string
        inputBinding:
            position: 4
    sample_name:
        type: string
        inputBinding:
            prefix: "-s"
outputs:
    annotated_expression_vcf:
        type: File
        outputBinding:
            glob: "annotated.expression.vcf.gz"

