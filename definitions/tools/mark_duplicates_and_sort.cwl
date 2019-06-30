#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Mark duplicates and Sort"

baseCommand: ["/bin/bash", "markduplicates_helper.sh"]
requirements:
    - class: ResourceRequirement
      coresMin: 8
      ramMin: 40000
    - class: DockerRequirement
      dockerPull: "mgibio/mark_duplicates-cwl:1.0.1"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'markduplicates_helper.sh'
        entry: |
            set -o pipefail
            set -o errexit

            declare MD_BARCODE_TAG
            if [ ! -z "$6" ]; then
              MD_BARCODE_TAG="BARCODE_TAG=$6"
            /usr/bin/java -Xmx16g -jar /opt/picard/picard.jar MarkDuplicates I=$1 O=/dev/stdout ASSUME_SORT_ORDER=$5 METRICS_FILE=$4 QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT "$MD_BARCODE_TAG" | /usr/bin/sambamba sort -t $2 -m 18G -o $3 /dev/stdin
            else
              /usr/bin/java -Xmx16g -jar /opt/picard/picard.jar MarkDuplicates I=$1 O=/dev/stdout ASSUME_SORT_ORDER=$5 METRICS_FILE=$4 QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT | /usr/bin/sambamba sort -t $2 -m 18G -o $3 /dev/stdin
            fi
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
    input_sort_order:
        type: string
        default: "queryname"
        inputBinding:
            position: 5
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
