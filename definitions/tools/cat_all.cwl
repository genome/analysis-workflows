#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/bin/cat']
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
arguments: [
    { shellQuote: false, valueFrom: "|" },
    "/bin/grep", "ChrID", "/dev/stdin"
]
stdout: "all_chromosome_pindel.head"
inputs:
    chromosome_pindel_outs:
        type: File[]
        inputBinding:
            position: -1 
outputs:
    all_chromosome_pindel_head:
        type: stdout
