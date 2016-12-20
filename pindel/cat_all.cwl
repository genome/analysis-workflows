#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/bin/cat']
requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
stdout: "all_interval_pindel.out"
inputs:
    interval_pindel_outs:
        type: File[]
        inputBinding:
            position: 1 
outputs:
    all_interval_pindel_out:
        type: stdout

