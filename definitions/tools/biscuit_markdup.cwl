#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool
label: "Biscuit dedup"
baseCommand: ["/bin/bash", "biscuit_markdup.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 24000
      coresMin: 4
    - class: DockerRequirement
      dockerPull: "mgibio/biscuit:0.3.8"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'biscuit_markdup.sh'
        entry: |
            set -eou pipefail

            cores=$1
            outdir="$2"
            bam="$3"

            /usr/bin/biscuit markdup "$bam" /dev/stdout | /usr/bin/sambamba sort -t $cores -m 15G -o "$outdir/markdup.bam" /dev/stdin
arguments:
    [{ valueFrom: $(runtime.cores), position: -3},
    { valueFrom: $(runtime.outdir), position: -2}
    ]
inputs:
    bam:
        type: File
        inputBinding:
            position: -1
outputs:
    markdup_bam:
        type: File
        outputBinding:
            glob: "markdup.bam"
