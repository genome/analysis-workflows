#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Xenosplit"

baseCommand: ["/bin/bash", "Xenosplit.sh"]

requirements:
    - class: ResourceRequirement
      ramMin: 20000
      tmpdirMin: 100000
    - class: DockerRequirement
      dockerPull: "mgibio/xenosplit:0.2"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'Xenosplit.sh'
        entry: |
            set -o pipefail
            set -o errexit
            
            python /opt/Xenosplit.py --pairedEnd --out graftOut.bam $1 $2
            python /opt/Xenosplit.py --count $1 $2 > goodnessOfMapping.txt

inputs:
    graftbam:
        type: File
        inputBinding:
            position: 1
    hostbam:
        type: File
        inputBinding:
            position: 2

outputs:
    graftOut:
        type: File
        outputBinding:
            glob: "graftOut.bam"
    goodnessOfMapping:
        type: File
        outputBinding:
            glob: "goodnessOfMapping.txt"
