#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'filter consensus reads'
baseCommand: ["/usr/local/bin/fgbio", "FilterConsensusReads"]
arguments:
    [ "--output", { valueFrom: "$(runtime.outdir)/consensus_filtered.bam"} ]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: quay.io/biocontainers/fgbio:1.3.0--0
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "--input"
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "--ref"
    min_reads:
        type: int[]
        inputBinding:
            prefix: "--min-reads"
    max_read_error_rate:
        type: float?
        inputBinding:
            prefix: "--max-read-error-rate"
    max_base_error_rate:
        type: float?
        inputBinding:
            prefix: "--max-base-error-rate"
    min_base_quality:
        type: int
        inputBinding:
            prefix: "--min-base-quality"
    max_no_call_fraction:
        type: float?
        inputBinding:
            prefix: "--max-no-call-fraction"
outputs:
    filtered_bam:
        type: File
        outputBinding:
            glob: "consensus_filtered.bam"
