#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/bin/sed', '/^\@/d']
requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
stdout: "interval_list.bed"
inputs:
    interval_list:
        type: File
        inputBinding:
            position: 1
outputs:
    interval_bed:
        type: stdout

