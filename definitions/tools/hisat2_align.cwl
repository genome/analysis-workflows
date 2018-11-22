#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "HISAT2: align"
baseCommand: ["/usr/bin/hisat2"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 16
    - class: DockerRequirement
      dockerPull: "mgibio/rnaseq"
arguments: [
    "-p", $(runtime.cores),
    "--dta",
    { shellQuote: false, valueFrom: "|" },
    "/usr/bin/sambamba", "view", "-S", "-f", "bam", "-l", "0", "/dev/stdin",
    { shellQuote: false, valueFrom: "|" },
    "/usr/bin/sambamba", "sort", "-t", $(runtime.cores), "-m", "8G", "-o", "$(runtime.outdir)/aligned.bam", "/dev/stdin"
]
inputs:
    reference_index:
        type: File
        secondaryFiles: [".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"]
        inputBinding:
            prefix: "-x"
            position: -3
    fastq1:
        type: File
        inputBinding:
            prefix: "-1"
            position: -2
    fastq2:
        type: File
        inputBinding:
            prefix: "-2"
            position: -1
    read_group_id:
        type: string
        inputBinding:
            prefix: "--rg-id"
            position: -5
    read_group_fields:
        type:
            type: array
            items: string
            inputBinding:
                prefix: "--rg"
        inputBinding:
            position: -4
    strand:
        type: string?
        inputBinding:
            valueFrom: |
                ${
                    if (inputs.strand) {
                        if (inputs.strand == 'first') {
                            return ['--rna-strandness RF'];
                        } else if (inputs.strand == 'second') {
                            return ['--rna-strandness FR'];
                        } else {
                            return [];
                        }
                    } else {
                            return []
                    }
                }
            position: -6
outputs:
    aligned_bam:
        type: File
        outputBinding:
            glob: "aligned.bam"
