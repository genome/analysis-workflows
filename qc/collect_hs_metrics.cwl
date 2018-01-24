#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect HS metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectHsMetrics"]
arguments:
    ["O=", { valueFrom: $(runtime.outdir)/HsMetrics.txt }]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
inputs:
    cram:
        type: File
        inputBinding:
            prefix: "I="
        secondaryFiles: [^.crai]
    reference:
        type: string
        inputBinding:
            prefix: "R="
    metric_accumulation_level:
        type: string
        inputBinding:
            prefix: "METRIC_ACCUMULATION_LEVEL="
    bait_intervals:
        type: File
        inputBinding:
            prefix: "BAIT_INTERVALS="
    target_intervals:
        type: File
        inputBinding:
            prefix: "TARGET_INTERVALS="
    per_target_coverage:
        type: boolean?
        inputBinding:
            prefix: "PER_TARGET_COVERAGE=PerTargetCoverage.txt"
    per_base_coverage:
        type: boolean?
        inputBinding:
            prefix: "PER_BASE_COVERAGE=PerBaseCoverage.txt"
outputs:
    hs_metrics:
        type: File
        outputBinding:
            glob: "HsMetrics.txt"
    per_target_coverage_metrics:
        type: File?
        outputBinding:
            glob: "PerTargetCoverage.txt"
    per_base_coverage_metrics:
        type: File?
        outputBinding:
            glob: "PerBaseCoverage.txt"
