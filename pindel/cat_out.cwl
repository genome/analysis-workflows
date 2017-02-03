#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/bin/cat']
requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
stdout: "per_interval_pindel.out"
inputs:
    pindel_outs:
        type: File[]
        inputBinding:
            position: 1 
outputs:
    pindel_out:
        type: stdout

