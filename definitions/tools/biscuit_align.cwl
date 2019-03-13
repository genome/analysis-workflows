#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Biscuit: align"
baseCommand: ["/bin/bash","biscuit_align.sh"]
requirements:
    - class: ResourceRequirement
      ramMin: 32000
      coresMin: 12
    - class: DockerRequirement
      dockerPull: "mgibio/biscuit:0.3.8"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'biscuit_align.sh'
        entry: |
            set -eou pipefail

            cores=$1
            outdir="$2"
            read_group_id="$3"
            reference_index="$4"
            fastq1="$5"
            fastq2="$6"

            /usr/bin/biscuit align -t $cores -M -R "$read_group_id" "$reference_index" "$fastq1" "$fastq2" | /usr/bin/sambamba view -S -f bam -l 0 /dev/stdin | /usr/bin/sambamba sort -t $cores -m 8G -o "$outdir/aligned.bam" /dev/stdin

arguments: [
    { valueFrom: $(runtime.cores), position: -9 },
    { valueFrom: $(runtime.outdir), position: -8 },
]
inputs:
    reference_index:
        type: string
        inputBinding:
            position: -3
    fastq1:
        type: File
        inputBinding:
            position: -2
    fastq2:
        type: File
        inputBinding:
            position: -1
    read_group_id:
        type: string
        inputBinding:
            position: -4
outputs:
    aligned_bam:
        type: File
        outputBinding:
            glob: "aligned.bam"
