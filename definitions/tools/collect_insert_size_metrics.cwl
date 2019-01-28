#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect insert size metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectInsertSizeMetrics"]
arguments:
    ["O=", { valueFrom: $(runtime.outdir)/$(inputs.cram.nameroot).InsertSizeMetrics.txt },
    "H=", { valueFrom: $(runtime.outdir)/$(inputs.cram.nameroot).InsertSizeHistogram.pdf }]
requirements:
    - class: ResourceRequirement
      ramMin: 18000
    - class: DockerRequirement
      dockerPull: "mgibio/picard-cwl:2.18.1"
inputs:
    cram:
        type: File
        inputBinding:
            prefix: "I="
        secondaryFiles: [^.crai]
    reference:
        type: string
        inputBinding:
            prefix: "REFERENCE_SEQUENCE="
    metric_accumulation_level:
        type: string
        inputBinding:
            prefix: "METRIC_ACCUMULATION_LEVEL="
outputs:
    insert_size_metrics:
        type: File
        outputBinding:
            glob: "$(inputs.cram.nameroot).InsertSizeMetrics.txt"
    insert_size_histogram:
        type: File
        outputBinding:
            glob: "$(inputs.cram.nameroot).InsertSizeHistogram.pdf"
