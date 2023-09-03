#! /usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: ["/bin/bash","directory_gatherer.sh"]

requirements:
    - class: DockerRequirement
      dockerPull: "ubuntu:focal"
    - class: ResourceRequirement
      ramMin: 1000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'directory_gatherer.sh'
        entry: |
            set -eou pipefail

            outdir="$1"
            files=("${@:2}")
            mkdir "$outdir"
            chmod -R 777 "$outdir"
            cp --recursive --preserve --no-clobber --target-directory "$outdir" "${files[@]}"

            exit 0

inputs:
    outdir:
        type: string
        inputBinding:
            position: 1
    directories:
         type: Directory[]
         inputBinding:
            position: 2
outputs:
    gathered_directory:
        type: Directory
        outputBinding:
            glob: "$(inputs.outdir)"

