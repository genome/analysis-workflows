#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect insert size metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectInsertSizeMetrics"]
arguments:
    ["O=", { valueFrom: $(runtime.outdir)/$(inputs.bam.nameroot).InsertSizeMetrics.txt },
    "H=", { valueFrom: $(runtime.outdir)/$(inputs.bam.nameroot).InsertSizeHistogram.pdf }]
requirements:
    - class: ResourceRequirement
      ramMin: 18000
    - class: DockerRequirement
      dockerPull: "mgibio/picard-cwl:2.18.1"
inputs:
    bam:
        type: File
        inputBinding:
            prefix: "I="
        secondaryFiles: [^.bai]
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
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
            glob: "$(inputs.bam.nameroot).InsertSizeMetrics.txt"
    insert_size_histogram:
        type: File
        outputBinding:
            glob: "$(inputs.bam.nameroot).InsertSizeHistogram.pdf"
