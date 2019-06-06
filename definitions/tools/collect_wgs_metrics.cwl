#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect WGS metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectWgsMetrics"]
arguments:
    ["O=", { valueFrom: $(runtime.outdir)/WgsMetrics.txt }]
requirements:
    - class: ResourceRequirement
      ramMin: 18000
    - class: DockerRequirement
      dockerPull: mgibio/picard-cwl:2.18.1
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "I="
        secondaryFiles: [^.bai]
    reference:
        type: string
        inputBinding:
            prefix: "R="
    intervals:
        type: File
        inputBinding:
            prefix: "INTERVALS="
outputs:
    wgs_metrics:
        type: File
        outputBinding:
            glob: "WgsMetrics.txt"
