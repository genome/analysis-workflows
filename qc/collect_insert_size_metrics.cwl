#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect insert size metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectInsertSizeMetrics"]
arguments:
    ["O=", { valueFrom: $(runtime.outdir)/InsertSizeMetrics.txt },
    "H=", { valueFrom: $(runtime.outdir)/InsertSizeHistogram.pdf }]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
inputs:
    cram:
        type: File
        inputBinding:
            prefix: "I="
        secondaryFiles: [^.crai]
outputs:
    insert_size_metrics:
        type: File
        outputBinding:
            glob: "InsertSizeMetrics.txt"
