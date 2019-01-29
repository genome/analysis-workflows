#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "run bam-readcount"

baseCommand: ["/usr/bin/python", "/usr/bin/bam_readcount_helper.py"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/bam_readcount_helper-cwl:1.0.0"
arguments: [
    { valueFrom: $(runtime.outdir), position: -3 }
]
stdout: $(inputs.sample)_bam_readcount.tsv
inputs:
    vcf:
        type: File
        inputBinding:
            position: -7
    sample:
        type: string
        inputBinding:
            position: -6
    reference_fasta:
        type: string
        inputBinding:
            position: -5
    bam:
        type: File
        inputBinding:
            position: -4
        secondaryFiles: [.bai]
    min_base_quality:
        type: int?
        default: 20
        inputBinding:
            position: -2
    min_mapping_quality:
        type: int?
        default: 0
        inputBinding:
            position: -1
outputs:
    snv_bam_readcount_tsv:
        type: File
        outputBinding:
            glob: "$(inputs.sample)_bam_readcount_snv.tsv"
    indel_bam_readcount_tsv:
        type: File
        outputBinding:
            glob: "$(inputs.sample)_bam_readcount_indel.tsv"
