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
            read_group_id="$2"
            reference_index="$3"
            fastq1="$4"
            if [[ $# -gt 5 ]];then
               echo "ERROR: too many arguments - were more than two fastq files provided to biscuit?"
               exit 1
            fi
            if [[ $# -gt 4 ]];then  #two fastqs
                fastq2=$5
                /usr/bin/biscuit align -t $cores -M -R "$read_group_id" "$reference_index" "$fastq1" "$fastq2" | /usr/bin/sambamba view -S -f bam -l 0 /dev/stdin | /usr/bin/sambamba sort -t $cores -m 8G -o "aligned.bam" /dev/stdin
            else #one fastq
                /usr/bin/biscuit align -t $cores -M -R "$read_group_id" "$reference_index" "$fastq1" | /usr/bin/sambamba view -S -f bam -l 0 /dev/stdin | /usr/bin/sambamba sort -t $cores -m 8G -o "aligned.bam" /dev/stdin
            fi

arguments: [
    { valueFrom: $(runtime.cores), position: 1 }
]
inputs:
    reference_index:
        type: string
        inputBinding:
            position: 3
    fastqs:
        type: File[]
        inputBinding:
            position: 4
    read_group_id:
        type: string
        inputBinding:
            position: 2
outputs:
    aligned_bam:
        type: File
        outputBinding:
            glob: "aligned.bam"
