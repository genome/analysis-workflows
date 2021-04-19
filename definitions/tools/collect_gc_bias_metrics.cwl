#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect gc bias metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectGcBiasMetrics"]
arguments:
    ["O=", { valueFrom: $(runtime.outdir)/$(inputs.sample_name).GcBiasMetrics.txt },
    "CHART=", { valueFrom: $(runtime.outdir)/$(inputs.sample_name).GcBiasMetricsChart.pdf },
    "S=", { valueFrom: $(runtime.outdir)/$(inputs.sample_name).GcBiasMetricsSummary.txt } ]
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
    metric_accumulation_level:
        type: string
        inputBinding:
            prefix: "METRIC_ACCUMULATION_LEVEL="
outputs:
    gc_bias_metrics:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).GcBiasMetrics.txt"
    gc_bias_metrics_chart:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).GcBiasMetricsChart.pdf"
    gc_bias_metrics_summary:
        type: File
        outputBinding:
            glob: "$(inputs.sample_name).GcBiasMetricsSummary.txt"
