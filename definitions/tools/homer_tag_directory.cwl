#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'Creating the homer tag directory'
baseCommand: ["/bin/bash", "/gscmnt/gc2708/info/medseq/apl/chipseq/src/bam2homer.sh"]
requirements:
    - class: DockerRequirement
      dockerPull: "chrisamiller/homer"
    - class: ResourceRequirement
      ramMin: 32000

inputs:
    tag_directory_name:
        type: string
        inputBinding:
            position: 1
    sam:
        type: File
        inputBinding:
            position: 2
    reference:
        type: string
        inputBinding:
            position: 3
outputs:
    tag_directory:
        type: Directory
        outputBinding:
            glob: "$(inputs.tag_directory_name)"
