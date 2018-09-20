#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/usr/bin/perl', '/usr/bin/intervals_to_bed.pl']
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/perl_helper-cwl:1.0.0"
stdout: "interval_list.bed"
inputs:
    interval_list:
        type: File
        inputBinding:
            position: 1
outputs:
    interval_bed:
        type: stdout

