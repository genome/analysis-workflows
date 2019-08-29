#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "bedGraph to bigwig conversion"
baseCommand: ["/bin/bash","bed2bigwig.sh"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 32000
      coresMin: 1
    - class: DockerRequirement
      dockerPull: "mgibio/bisulfite:v1.4"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'bed2bigwig.sh'
        entry: |
            set -eou pipefail
                
                for FILE in ${@:2}
                do
                    filename=`basename $FILE .bedgraph`
                    /usr/bin/bedGraphToBigWig $FILE $1 $filename.bw
                done
inputs:
    methylation_bedgraph:
        type: File[]
        inputBinding:
            position: 2
    reference_sizes:
        type: File
        inputBinding:
            position: 1
outputs:
    methylation_bigwig:
        type: File[]
        outputBinding:
            glob: "*.bw"
