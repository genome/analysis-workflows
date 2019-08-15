#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'Creating the homer tag directory'
baseCommand: ["/bin/bash", "homer_tag_directory.sh"]
requirements:
    - class: DockerRequirement
      dockerPull: "chrisamiller/homer"
    - class: ResourceRequirement
      ramMin: 32000
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'homer_tag_directory.sh'
        entry: |
                set -o pipefail
                set -o errexit

                name=$1
                bam=$2
                reference_fasta=$3  #assumes that there's a .genome file next to the reference .fa
                export PATH=$PATH:/opt/homer/bin
                if [[ ! -d $name/homer ]];then
                    mkdir -p $name/homer
                fi
                echo "creating tagDir"
                makeTagDirectory $name/homer $bam

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
