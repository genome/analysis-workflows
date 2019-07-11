#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "align with bwa_mem and tag"

baseCommand: ["/bin/bash", "fastq_alignment_helper.sh"]
requirements:
    - class: ResourceRequirement
      coresMin: 8
      ramMin: 20000
    - class: DockerRequirement
      dockerPull: "mgibio/alignment_helper-cwl:1.0.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'fastq_alignment_helper.sh'
        entry: |
            set -o pipefail
            set -o errexit
            set -o nounset

            FASTQ1="$1"
            FASTQ2="$2"
            READGROUP="$3"
            REFERENCE="$4"
            NTHREADS="$5"

            /usr/local/bin/bwa mem -K 100000000 -t "$NTHREADS" -Y -p -R "$READGROUP" "$REFERENCE" "$FASTQ1" "$FASTQ2" | /usr/local/bin/samblaster -a --addMateTags | /opt/samtools/bin/samtools view -b -S /dev/stdin
stdout: "refAlign.bam"
arguments:
    - position: 5
      valueFrom: $(runtime.cores)
inputs:
    fastq_1:
        type: File
        inputBinding:
            position: 1
        doc: 'fastq of sequence data (first read if paired end)'
    fastq_2:
        type: File
        inputBinding:
            position: 2
        doc: 'second read fastq for paired end data'
    readgroup:
        type: string
        inputBinding:
            position: 3
        doc: 'readgroup information for the BAM header'
    reference:
        type: string
        inputBinding:
            position: 4
        doc: 'bwa-indexed reference file'
outputs:
    aligned_bam:
        type: stdout
