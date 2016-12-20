#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/bin/cat']
requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
stdout: "all_interval_pindel.head"
inputs:
    all_interval_pindel_heads:
        type: File[]
        inputBinding:
            position: 1 
outputs:
    all_interval_pindel_head:
        type: stdout

