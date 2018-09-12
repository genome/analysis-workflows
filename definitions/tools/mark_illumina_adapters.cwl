#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'fastq to bam'
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/picard/picard.jar", "MarkIlluminaAdapters"]
arguments:
    ["OUTPUT=", { valueFrom: $(runtime.outdir)/marked.bam },
    "METRICS=", { valueFrom: $(runtime.outdir)/adapter_metrics.txt }]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "mgibio/dna-alignment"
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "INPUT="
            position: 1
outputs:
    marked_bam:
        type: File
        outputBinding:
            glob: "marked.bam"
    metrics:
        type: File
        outputBinding:
            glob: "adapter_metrics.txt"
