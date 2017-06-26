#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "run bam-readcount"

baseCommand: ["/usr/bin/python", "/usr/bin/bam_readcount_helper.py"]
arguments:
    - position: 5
      valueFrom: $(runtime.outdir)
inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
    sample:
        type: string
        inputBinding:
            position: 2
    reference_fasta:
        type: string
        inputBinding:
            position: 3
    bam:
        type: File
        inputBinding:
            position: 4
        secondaryFiles: [.bai]
outputs:
    snv_bam_readcount:
        type: File?
        outputBinding:
            glob: bam_readcount_snv.tsv
    indel_bam_readcount:
        type: File?
        outputBinding:
            glob: bam_readcount_indel.tsv
