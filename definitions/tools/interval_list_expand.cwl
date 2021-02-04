#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "expand interval list regions by a given number of basepairs"

baseCommand: ["/usr/bin/java", "-Xmx3g", "-jar", "/usr/picard/picard.jar", "IntervalListTools"]

arguments:
    [valueFrom: "OUTPUT=$(runtime.outdir)/$(inputs.interval_list.nameroot).expanded.interval_list", "UNIQUE=TRUE"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "broadinstitute/picard:2.23.6"
inputs:
    interval_list:
        type: File
        inputBinding:
            prefix: "INPUT="
            separate: false
    roi_padding:
        type: int
        inputBinding:
            prefix: "PADDING="
            separate: false
outputs:
    expanded_interval_list:
        type: File
        outputBinding:
            glob: "$(inputs.interval_list.nameroot).expanded.interval_list"
