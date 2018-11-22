#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "StringTie"
baseCommand: ["/usr/bin/stringtie"]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 12
    - class: DockerRequirement
      dockerPull: "mgibio/rnaseq"
arguments: [
    "-o", "$(runtime.outdir)/stringtie_transcripts.gtf",
    "-A", "$(runtime.outdir)/stringtie_gene_expression.tsv",
    "-p", $(runtime.cores),
    "-e"
]
inputs:
    strand:
        type: string?
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
