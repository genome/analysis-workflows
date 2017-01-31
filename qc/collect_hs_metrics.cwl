#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect HS metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectHsMetrics"]
arguments:
    ["O=", { valueFrom: $(runtime.outdir)/HsMetrics.txt }]
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
    bait_intervals:
        type: File
        inputBinding:
            prefix: "BAIT_INTERVALS="
    target_intervals:
        type: File
        inputBinding:
            prefix: "TARGET_INTERVALS="
outputs:
    hs_metrics:
        type: File
        outputBinding:
            glob: "HsMetrics.txt"
