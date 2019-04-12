#! /usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["/bin/bash","directory_gatherer.sh"]

requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "ubuntu:xenial"
    - class: ResourceRequirement
      ramMin: 8000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'directory_gatherer.sh'
        entry: |
            outdir="$1"
            files="${@:-1}"
            cp $files $outdir

inputs:
    outdir:
        type: string
        inputBinding:
            position: 1
    files:
        type: File[]
        inputBinding:
            position: 2
outputs:
    gathered_directory:
        type: Directory
        outputBinding:
            glob: "$(inputs.outdir)"
