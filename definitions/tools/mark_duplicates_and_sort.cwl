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
      dockerPull: "mgibio/mark_duplicates-cwl:2.0.0"
    - class: InitialWorkDirRequirement
      listing:
      - entryname: 'markduplicates_helper.sh'
        entry: |
            set -o pipefail
            set -o errexit

            /usr/bin/java -Xmx16g -jar /opt/picard/picard.jar MarkDuplicates I=$1 O=/dev/stdout ASSUME_SORT_ORDER=$5 METRICS_FILE=$4 QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT REFERENCE_SEQUENCE=$6 | /opt/samtools/bin/samtools sort -@ $2 -m 4G --reference "$6" -o "$3" -O cram /dev/stdin

arguments:
    - position: 2
      valueFrom: "$(runtime.cores)"
    - position: 4
      valueFrom: "$(inputs.cram.nameroot).mark_dups_metrics.txt"
inputs:
    cram:
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
        default: 'MarkedSorted'
        inputBinding:
            position: 3
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai]
        inputBinding:
            position: 6

outputs:
    sorted_cram:
        type: File
        outputBinding:
            glob: $(inputs.output_name).cram
        secondaryFiles: [.crai]
    metrics_file:
        type: File
        outputBinding:
            glob: "$(inputs.cram.nameroot).mark_dups_metrics.txt"
