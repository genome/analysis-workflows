#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "downsample alignment"
baseCommand: ["/usr/bin/java", "-Xmx16g", "-jar", "/usr/picard/picard.jar", "DownsampleSam"]
arguments:
    ["OUTPUT=", { valueFrom: $(runtime.outdir)/downsampled.cram }]
requirements:
    - class: ResourceRequirement
      ramMin: 16000
inputs:
    cram:
        type: File
        inputBinding:
            prefix: "INPUT="
        secondaryFiles: [^.crai]
    probability:
        type: float
        inputBinding:
            prefix: "PROBABILITY="
outputs:
    downsampled_cram:
        type: File
        outputBinding:
            glob: "downsampled.cram"
