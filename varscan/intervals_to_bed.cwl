#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
baseCommand: ['/bin/bash', 'helper.sh']
requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
    - class: InitialWorkDirRequirement
      listing:
          - entryname: 'helper.sh'
            entry: |
                /bin/grep -v '^@' $1 | awk '{print $1 "\\t" $2-1 "\\t" $3}'
stdout: "interval_list.bed"
inputs:
    interval_list:
        type: File
        inputBinding:
            position: 1
outputs:
    interval_bed:
        type: stdout

