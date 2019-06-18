#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/bin/cat']
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
    - class: ResourceRequirement
      ramMin: 4000
arguments: [
    { shellQuote: false, valueFrom: "|" },
    "/bin/grep", "ChrID", "/dev/stdin"
]
stdout: "all_region_pindel.head"
inputs:
    region_pindel_outs:
        type: File[]
        inputBinding:
            position: -1 
outputs:
    all_region_pindel_head:
        type: stdout
