#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/usr/bin/perl', '/usr/bin/split_interval_list_helper.pl']
requirements:
    - class: ResourceRequirement
      ramMin: 6000
    - class: DockerRequirement
      dockerPull: mgibio/cle
arguments:
    [{ valueFrom: OUTPUT=$(runtime.outdir) }]
inputs:
    interval_list:
        type: File
        inputBinding:
            prefix: "INPUT="
            separate: false
            position: 1
    scatter_count:
        type: int
        inputBinding:
            prefix: "SCATTER_COUNT="
            separate: false
            position: 2
outputs:
    split_interval_lists:
        type: File[]
        outputBinding:
            glob: "*.interval_list"
