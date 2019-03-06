#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Biscuit: align"
baseCommand: ["/usr/bin/biscuit","align"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 32000
      coresMin: 12
    - class: DockerRequirement
      dockerPull: "mgibio/bisulfite:v1.3"
arguments: [
    { valueFrom: "-t", position: -10 },
    { valueFrom: $(runtime.cores), position: -9 },
    { valueFrom: "-M", position: -8 },
    { shellQuote: false, valueFrom: "|" },
    "/usr/bin/sambamba", "view", "-S", "-f", "bam", "-l", "0", "/dev/stdin",
    { shellQuote: false, valueFrom: "|" },
    "/usr/bin/sambamba", "sort", "-t", $(runtime.cores), "-m", "8G", "-o", "$(runtime.outdir)/aligned.bam", "/dev/stdin"
]
inputs:
    reference_index:
        type: string
        inputBinding:
            position: -3
    fastq1:
        type: File
        inputBinding:
            position: -2
    fastq2:
        type: File
        inputBinding:
            position: -1
    read_group_id:
        type: string
        inputBinding:
            prefix: "-R"
            position: -4
outputs:
    aligned_bam:
        type: File
        outputBinding:
            glob: "aligned.bam"
