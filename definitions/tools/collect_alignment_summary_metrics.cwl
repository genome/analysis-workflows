#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "collect alignment summary metrics"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "CollectAlignmentSummaryMetrics"]
arguments:
    ["OUTPUT=", { valueFrom: $(runtime.outdir)/$(inputs.bam.nameroot).AlignmentSummaryMetrics.txt }]
requirements:
    - class: ResourceRequirement
      ramMin: 18000
    - class: DockerRequirement
      dockerPull: "broadinstitute/picard:2.23.6"
    - class: InitialWorkDirRequirement
      listing:
      - $(inputs.ref_fai)
      - $(inputs.ref_dict)
      - $(inputs.reference)

inputs:
    ref_fai:
        type: File
    ref_dict:
        type: File

    bam:
        type: File
        inputBinding:
            prefix: "INPUT="
        secondaryFiles: [^.bai]
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            valueFrom: $(self.basename)
            prefix: "REFERENCE_SEQUENCE="
    metric_accumulation_level:
        type: string
        inputBinding:
            prefix: "METRIC_ACCUMULATION_LEVEL="
outputs:
    alignment_summary_metrics:
        type: File
        outputBinding:
            glob: "$(inputs.bam.nameroot).AlignmentSummaryMetrics.txt"
