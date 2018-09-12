#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: [
    "bash", "-c", "/bin/grep -v '^@' $1 | /usr/bin/cut -f 1 | /usr/bin/sort | /usr/bin/uniq", "--"
]
requirements:
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
stdout: "chromosome.list"
inputs:
    interval_list:
        type: File
        inputBinding:
            position: 1
outputs:
    chromosome_list:
        type: 
            type: array
            items: string
        outputBinding:
            glob: "chromosome.list"
            loadContents: true
            outputEval: $( self[0].contents.split("\n").slice(0, -1))


