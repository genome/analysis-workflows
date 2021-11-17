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
            files="${@:2}"
            mkdir $outdir
            chmod -R 777 $outdir
            cp --recursive --preserve --no-clobber --target-directory "$outdir" "$files"

            exit 0

inputs:
    outdir:
        type: string
        inputBinding:
            position: 1
    files:
        type: File[]
        inputBinding:
            position: 2
            valueFrom: |
              ${
                var results = []
                for(var i=0; i<self.length; i++){
                  results.push(self[i])
                  if(self[i].hasOwnProperty('secondaryFiles')){
                    for(var j=0; j<self[i].secondaryFiles.length; j++){
                      results.push(self[i].secondaryFiles[j])
                    }
                  }
                }
                return results
              }
    directory:
         type: Directory?
         inputBinding:
            position: 3
outputs:
    gathered_directory:
        type: Directory
        outputBinding:
            glob: "$(inputs.outdir)"

