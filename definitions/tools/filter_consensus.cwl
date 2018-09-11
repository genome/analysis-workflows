#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'filter consensus reads'
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/fgbio-0.5.0.jar", "FilterConsensusReads"]
arguments:
    ["--min-reads", "10", "5", "3", "--max-read-error-rate", "0.05",
    "--max-base-error-rate", "0.1", "--min-base-quality", "50",
    "--max-no-call-fraction", "0.05",
    "--output", { valueFrom: "$(runtime.outdir)/consensus_filtered.bam"} ]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: mgibio/dna-alignment
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "--input"
    reference:
        type: string
        inputBinding:
            prefix: "--ref"
outputs:
    filtered_bam:
        type: File
        outputBinding:
            glob: "consensus_filtered.bam"
