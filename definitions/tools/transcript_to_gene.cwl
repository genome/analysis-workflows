#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Kallisto: TranscriptToGene"
baseCommand: ["/usr/local/bin/Rscript"]
arguments: [
    "/usr/src/transcript_to_gene.R"
]
requirements:
    - class: ResourceRequirement
      ramMin: 2000
      coresMin: 1
    - class: DockerRequirement
      dockerPull: "mgibio/rnaseq"
inputs:
    gene_transcript_lookup_table:
        type: File
        inputBinding:
            position: 1
    transcript_table_h5:
        type: File
        inputBinding:
            position: 2
outputs:
    gene_abundance:
        type: File
        outputBinding:
            glob: "gene_abundance.tsv"
