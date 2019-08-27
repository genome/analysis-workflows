#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Staged Renamer"
doc: "Renames a file by staging and then `mv`ing it.  A workaround for workflow engines that don't support rename.cwl.  If running in cwltool, use the other one instead."
requirements:
    - class: ResourceRequirement
      ramMin: 4000
      coresMin: 1
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
