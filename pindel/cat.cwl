#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/bin/cat']
requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
stdout: "pindel.out"
inputs:
    deletion_out:
        type: File
        inputBinding:
            position: 1
    insertion_out:
        type: File
        inputBinding:
            position: 2 
outputs:
    pindel_out:
        type: stdout

