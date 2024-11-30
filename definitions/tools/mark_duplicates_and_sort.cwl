#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Mark duplicates and Sort"

baseCommand: ["/bin/bash", "markduplicates_helper.sh"]
requirements:
    - class: ResourceRequirement
      coresMin: 16
      ramMin: 40000
    - class: DockerRequirement
      dockerPull: "quay.io/biocontainers/sambamba:0.8.2--h98b6b92_2"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'markduplicates_helper.sh'
        entry: |
            set -o pipefail
            set -o errexit

            CORES="$2"
            CORES_PER_JOB=`perl -E 'my $x = int($ARGV[0]/2); say($x < 1? 1 : $x)'` $CORES

            sambamba markdup -l 0 -t $CORES_PER_JOB "$1" /dev/stdout 2> "$4" \
              | sambamba sort -t $CORES_PER_JOB -m 16G -o "$3" /dev/stdin
arguments:
    - position: 2
      valueFrom: "$(runtime.cores)"
    - position: 4
      valueFrom: "$(inputs.bam.nameroot).mark_dups_metrics.txt"
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
    output_name:
        type: string?
        default: 'MarkedSorted.bam'
        inputBinding:
            position: 3

outputs:
    sorted_bam:
        type: File
        outputBinding:
            glob: $(inputs.output_name)
        secondaryFiles: [.bai]
    metrics_file:
        type: File
        outputBinding:
            glob: "$(inputs.bam.nameroot).mark_dups_metrics.txt"
