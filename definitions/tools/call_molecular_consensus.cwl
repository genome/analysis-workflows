#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'call molecular consensus'
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/fgbio-0.5.0.jar", "CallMolecularConsensusReads"]
arguments:
    ["--error-rate-pre-umi", "45", "--error-rate-post-umi", "30", "--min-input-base-quality", "30", "--min-reads", "1",
    "--output", { valueFrom: "$(runtime.outdir)/consensus_unaligned.bam"} ]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: mgibio/dna-alignment:1.0.0
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "--input"
outputs:
    consensus_bam:
        type: File
        outputBinding:
            glob: "consensus_unaligned.bam"
