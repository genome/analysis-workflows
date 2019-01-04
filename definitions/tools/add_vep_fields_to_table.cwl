#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "add VEP annotation to report"

baseCommand: ["/usr/bin/python", "/usr/bin/add_annotations_to_table_helper.py"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/annotation_table-cwl:1.0.3"
arguments:
    - position: 4
      valueFrom: $(runtime.outdir)
inputs:
    tsv:
        type: File
        inputBinding:
            position: 1
    vcf:
        type: File
        inputBinding:
            position: 2
    vep_fields:
        type:
            type: array
            items: string
        inputBinding:
            prefix: ""
            itemSeparator: ","
            separate: false
            position: 3
outputs:
    annotated_variants_tsv:
        type: File
        outputBinding:
            glob: variants.annotated.tsv
