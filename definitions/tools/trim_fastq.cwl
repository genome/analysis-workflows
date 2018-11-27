#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Trim FASTQ (flexbar)"
baseCommand: ['/opt/flexbar/flexbar']
arguments: [
    "--target", {valueFrom: "$(runtime.outdir)/trimmed_read"},
    "--threads", {valueFrom: "$(runtime.cores)"}
]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
      tmpdirMin: 25000
      coresMin: 4
    - class: DockerRequirement
      dockerPull: "mgibio/bisulfite"
inputs:
    adapters:
        type: File
        inputBinding:
            prefix: "--adapters"
            position: 1
    adapter_trim_end:
        type: string
        inputBinding:
            prefix: "--adapter-trim-end"
            position: 2
    adapter_min_overlap:
        type: int
        inputBinding:
            prefix: "--adapter-min-overlap"
            position: 3
    max_uncalled:
        type: int
        inputBinding:
            prefix: "--max-uncalled"
            position: 5
    min_readlength:
        type: int
        inputBinding:
            prefix: "--min-read-length"
            position: 6
    reads1:
        type: File
        inputBinding:
            prefix: "--reads"
            position: 7
    reads2:
        type: File
        inputBinding:
            prefix: "--reads2"
            position: 8
outputs:
    fastq1:
        type: File
        outputBinding:
            glob: "trimmed_read_1.fastq"
    fastq2:
        type: File
        outputBinding:
            glob: "trimmed_read_2.fastq"
    fastqs:
        type: File[]
        outputBinding:
            glob: "trimmed_read_*.fastq"
