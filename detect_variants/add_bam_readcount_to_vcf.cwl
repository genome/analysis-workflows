#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "add bam_readcount info to vcf"

baseCommand: ["/usr/bin/python3", "/usr/bin/add_bam_readcount_to_vcf_helper.py"]
arguments:
    [$(runtime.outdir)]
inputs:
    vcf:
        type: File
        inputBinding:
            position: -3
    bam_readcount_tsvs:
        type:
            type: array
            items: File
        inputBinding:
            prefix: ""
            itemSeparator: ","
            separate: false
            position: -2
    sample_names:
        type:
            type: array
            items: string
        inputBinding:
            prefix: ""
            itemSeparator: ","
            separate: false
            position: -1
outputs:
    annotated_bam_readcount_vcf:
        type: File
        outputBinding:
            glob: "annotated.bam_readcount.vcf.gz"

