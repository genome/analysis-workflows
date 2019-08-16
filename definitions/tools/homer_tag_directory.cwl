#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'Creating the homer tag directory'
baseCommand: ["/bin/bash", "homer_tag_directory.sh"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/homer:0.1"
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
                reference_fasta=$2

                export PATH=$PATH:/opt/homer/bin
                if [[ ! -d $name/homer ]];then
                    mkdir -p $name/homer
                fi
                echo "creating tagDir"
                makeTagDirectory $name/homer $bam

inputs:
    sam:
        type: File
        inputBinding:
            position: 1
    reference:
        type: string
        inputBinding:
            position: 2
outputs:
    tag_directory:
        type: Directory
        outputBinding:
            glob: "homer_tag_directory"
