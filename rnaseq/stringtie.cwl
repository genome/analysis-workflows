#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "StringTie"
baseCommand: ["/usr/bin/stringtie"]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      coresMin: 16
arguments: [
    "-o", "$(runtime.outdir)/stringtie_transcripts.gtf",
    "-A", "$(runtime.outdir)/stringtie_gene_expression.tsv",
    "-n", $(runtime.cores)
]
inputs:
    reference_annotation:
        type: File
        inputBinding:
            prefix: "-G"
            position: 1
    sample_name:
        type: string
        inputBinding:
            prefix: "-l"
            position: 2
    bam:
        type: File
        inputBinding:
            position: 3
outputs:
    gtf:
        type: File
        outputBinding:
            glob: stringtie_transcripts.gtf
    tsv:
        type: File
        outputBinding:
            glob: stringtie_gene_expression.tsv
