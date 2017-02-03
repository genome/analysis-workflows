#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/bin/grep', 'ChrID']
requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
stdout: "all_interval_pindel.head"
inputs:
    pindel_output:
        type: File
        inputBinding:
            position: 1
outputs:
    pindel_head:
        type: stdout

