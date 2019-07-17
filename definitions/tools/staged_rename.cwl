#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "A crude file-renamer to work around workflow engines that don't support rename.cwl"
requirements:
    - class: DockerRequirement
      dockerPull: ubuntu:bionic
    - class: InitialWorkDirRequirement
      listing:
          - $(inputs.original)
baseCommand: ["/bin/mv"]
arguments: [$(inputs.original.basename), $(inputs.name)]
inputs:
    original:
        type: File
    name:
        type: string
outputs:
    replacement:
        type: File
        outputBinding:
            glob: $(inputs.name)
