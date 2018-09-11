#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'group reads by umi'
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/fgbio-0.5.0.jar", "GroupReadsByUmi"]
arguments:
    ["--strategy", "paired", "--assign-tag", "MI", "--raw-tag", "RX", "--min-map-q", "10", "--edits", "1",
    "--output", { valueFrom: "$(runtime.outdir)/umi_grouped.bam"} ]
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
outputs:
    grouped_bam:
        type: File
        outputBinding:
            glob: "umi_grouped.bam"
