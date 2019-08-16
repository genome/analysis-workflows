#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'Creating the homer tag directory'
baseCommand: ["/bin/bash", "homer_tag_directory.sh"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/homer:4.10"
    - class: ResourceRequirement
      ramMin: 32000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'homer_tag_directory.sh'
        entry: |
                set -o pipefail
                set -o errexit

                name="homer_tag_directory"
                bam=$1
                mkdir -p $name
                echo "creating tagDir"
                /opt/homer/bin/makeTagDirectory $name $bam

inputs:
    sam:
        type: File
        inputBinding:
            position: 1
outputs:
    tag_directory:
        type: Directory
        outputBinding:
            glob: "homer_tag_directory"
