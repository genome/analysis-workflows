#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "mark duplicates and sort"

baseCommand: ["/bin/bash", "markduplicates_helper.sh"]
requirements:
    - class: ResourceRequirement
      coresMin: 2
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: "mgibio/mark_duplicates-cwl:1.0.1"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'markduplicates_helper.sh'
        entry: |
            set -o pipefail
            set -o errexit

            declare MD_BARCODE_TAG
            if [ ! -z "${5}" ]; then
              MD_BARCODE_TAG="BARCODE_TAG=${5}"
            /usr/bin/java -Xmx16g -jar /opt/picard/picard.jar MarkDuplicates I=$1 O=/dev/stdout ASSUME_SORT_ORDER=$5 METRICS_FILE=$4 QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT "${MD_BARCODE_TAG}" | /usr/bin/sambamba sort -t $2 -m 18G -o $3 /dev/stdin
            else
              /usr/bin/java -Xmx16g -jar /opt/picard/picard.jar MarkDuplicates I=$1 O=/dev/stdout ASSUME_SORT_ORDER=$5 METRICS_FILE=$4 QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT | /usr/bin/sambamba sort -t $2 -m 18G -o $3 /dev/stdin
            fi
arguments:
    - position: 2
      valueFrom: $(runtime.cores)
    - position: 3
      valueFrom: $(runtime.outdir)/MarkedSorted.bam
    - position: 4
      valueFrom: "$(inputs.bam.nameroot).mark_dups_metrics.txt"
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
    markdup_assume_sort_order:
        type: string
        inputBinding:
            position: 5
outputs:
    sorted_bam:
        type: File
        outputBinding:
            glob: "MarkedSorted.bam"
        secondaryFiles: [.bai]
    metrics_file:
        type: File
        outputBinding:
            glob: "$(inputs.bam.nameroot).mark_dups_metrics.txt"
