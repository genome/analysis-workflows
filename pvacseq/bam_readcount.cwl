#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "run bam-readcount"

baseCommand: ["/usr/bin/python", "/usr/bin/bam_readcount_helper.py"]
requirements:
    - class: ShellCommandRequirement
arguments: [
    $(runtime.outdir),
    { valueFrom: " && ", shellQuote: false },
    "/bin/cat", "$(runtime.outdir)/$(inputs.sample)_bam_readcount_snv.tsv", "$(runtime.outdir)/$(inputs.sample)_bam_readcount_indel.tsv"
]
stdout: $(inputs.sample)_bam_readcount.tsv
inputs:
    vcf:
        type: File
        inputBinding:
            position: -6
    sample:
        type: string
        inputBinding:
            position: -5
    reference_fasta:
        type: string
        inputBinding:
            position: -4
    bam:
        type: File
        inputBinding:
            position: -3
        secondaryFiles: [.bai]
    min_mapping_quality:
        type: int?
        default: 0
        inputBinding:
            position: -2
    min_base_quality:
        type: int?
        default: 20
        inputBinding:
            position: -1
outputs:
    bam_readcount_tsv:
        type: stdout
