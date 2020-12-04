#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect WGS metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectWgsMetrics"]
arguments:
    ["O=", { valueFrom: $(runtime.outdir)/$(inputs.sample_name).WgsMetrics.txt }]
requirements:
    - class: ResourceRequirement
      ramMin: 18000
    - class: DockerRequirement
      dockerPull: broadinstitute/picard:2.23.6
inputs:
    sample_name:
        type: string?
        default: "final"
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
            prefix: "R="
    intervals:
        type: File?
        inputBinding:
            prefix: "INTERVALS="
    minimum_mapping_quality:
        type: int?
        doc: 'Minimum mapping quality for a read to contribute coverage.'
        inputBinding:
            prefix: "MINIMUM_MAPPING_QUALITY="
    minimum_base_quality:
        type: int?
        doc: 'Minimum base quality for a base to contribute coverage. (N bases are considered to have negative infinite quality and will be excluded regardless of this value.)'
        inputBinding:
            prefix: "MINIMUM_BASE_QUALITY="
outputs:
    wgs_metrics:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).WgsMetrics.txt"
