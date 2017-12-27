#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Kallisto: TranscriptToGene"
baseCommand: ["/usr/local/bin/Rscript"]
arguments: [
    "/usr/src/transcript_to_gene.R"
]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 2000
      coresMin: 16
inputs:
    gene_transcript_lookup_table:
        type: File
        inputBinding:
            position: 1
    transcript_table:
        type: File
        inputBinding:
            position: 2
outputs:
    gene_abundances:
        type: File
        outputBinding:
            glob: "gene_abundances.tsv"
    gene_counts:
        type: File
        outputBinding:
            glob: "gene_counts.tsv"
    gene_lengths:
        type: File
        outputBinding:
            glob: "gene_lengths.tsv"
