#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/usr/bin/perl', '/usr/bin/interval_to_bed_split.pl']
requirements:
    - class: ResourceRequirement
      ramMin: 6000
    - class: DockerRequirement
      dockerPull: mgibio/perl_helper-cwl:1.1.0
arguments: [$(runtime.outdir)]
inputs:
    interval_list:
        type: File
        inputBinding:
            position: 1
    scatter_count:
        type: int
        inputBinding:
            position: 2
outputs:
    split_beds:
        type: File[]
        outputBinding:
            glob: "*.interval.bed"
