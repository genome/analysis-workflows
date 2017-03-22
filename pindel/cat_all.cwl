#! /usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/bin/cat']
stdout: "all_chromosome_pindel.out"
inputs:
    chromosome_pindel_outs:
        type: File[]
        inputBinding:
            position: 1 
outputs:
    all_chromosome_pindel_out:
        type: stdout

