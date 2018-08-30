#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect HS metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectHsMetrics"]
arguments:
    ["O=", { valueFrom: $(runtime.outdir)/$(inputs.output_prefix)-HsMetrics.txt }]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
    - class: InlineJavascriptRequirement
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
            prefix: "PER_TARGET_COVERAGE="
            valueFrom: |
                        ${
                            if(self) {
                                return inputs.output_prefix + "-PerTargetCoverage.txt"
                            } else {
                                return false;
                            }
                        }
    per_base_coverage:
        type: boolean?
        inputBinding:
            prefix: "PER_BASE_COVERAGE="
            valueFrom: |
                        ${
                            if(self) {
                                return inputs.output_prefix + "-PerBaseCoverage.txt"
                            } else {
                                return false;
                            }
                        }
    output_prefix:
        type: string?
        default: "out"
    minimum_mapping_quality:
        type: int?
        default: 1
        inputBinding:
            prefix: "MINIMUM_MAPPING_QUALITY="
    minimum_base_quality:
        type: int?
        default: 1
        inputBinding:
            prefix: "MINIMUM_BASE_QUALITY="
outputs:
    hs_metrics:
        type: File
        outputBinding:
            glob: "$(input.output_prefix)-HsMetrics.txt"
    per_target_coverage_metrics:
        type: File?
        outputBinding:
            glob: "$(input.output_prefix)-PerTargetCoverage.txt"
    per_base_coverage_metrics:
        type: File?
        outputBinding:
            glob: "$(input.output_prefix)-PerBaseCoverage.txt"
