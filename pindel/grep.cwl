#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/bin/grep', 'ChrID']
stdout: "all_chromosome_pindel.head"
inputs:
    pindel_output:
        type: File
        inputBinding:
            position: 1
outputs:
    pindel_head:
        type: stdout

