#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "StringTie"
baseCommand: ["/usr/local/bin/stringtie"]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 12
    - class: DockerRequirement
      dockerPull: "quay.io/biocontainers/stringtie:2.1.4--h7e0af3c_0"
    - class: StepInputExpressionRequirement

arguments: [
    "-o", "$(runtime.outdir)/stringtie_transcripts.gtf",
    "-A", "$(runtime.outdir)/stringtie_gene_expression.tsv",
    "-p", $(runtime.cores),
    "-e"
]
inputs:
    strand:
        type:
          - "null"
          - type: enum
            symbols: ["first", "second", "unstranded"]
        inputBinding:
            valueFrom: |
                ${
                    if (inputs.strand) {
                        if (inputs.strand == 'first') {
                            return ['--rf'];
                        } else if (inputs.strand == 'second') {
                            return ['--fr'];
                        } else {
                            return [];
                        }
                    } else {
                            return []
                    }
                }
            position: 1
    reference_annotation:
        type: File
        inputBinding:
            prefix: "-G"
            position: 2
    sample_name:
        type: string
        inputBinding:
            prefix: "-l"
            position: 3
    bam:
        type: File
        inputBinding:
            position: 4
outputs:
    transcript_gtf:
        type: File
        outputBinding:
            glob: stringtie_transcripts.gtf
    gene_expression_tsv:
        type: File
        outputBinding:
            glob: stringtie_gene_expression.tsv
