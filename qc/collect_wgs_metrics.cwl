#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect WGS metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectWgsMetrics"]
arguments:
    ["O=", { valueFrom: $(runtime.outdir)/WgsMetrics.txt }]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/picard-2.4.1-r:2"
    - class: ResourceRequirement
      ramMin: 16000
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "I="
        secondaryFiles: [^.bai]
    reference:
        type: File
        inputBinding:
            prefix: "R="
        secondaryFiles: [.fai]
    intervals:
        type: File
        inputBinding:
            prefix: "INTERVALS="
outputs:
    wgs_metrics:
        type: File
        outputBinding:
            glob: "WgsMetrics.txt"
