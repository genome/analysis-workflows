#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "mark duplicates and sort"

baseCommand: ["/bin/bash", "helper.sh"]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/sort-mark-duplicates:2"
    - class: ResourceRequirement
      coresMin: 8
      ramMin: 18000
    - class: InitialWorkDirRequirement
      listing:
          - entryname: 'helper.sh'
            entry: |
                /usr/bin/java -Xmx16g -jar /usr/picard/picard.jar MarkDuplicates I=$1 O=/dev/stdout ASSUME_SORT_ORDER=queryname METRICS_FILE=mark_dups_metrics.txt QUIET=true COMPRESSION_LEVEL=0 VALIDATION_STRINGENCY=LENIENT | /usr/local/bin/sambamba sort -t $2 -m 18G -o $3 /dev/stdin
arguments:
    - position: 2
      valueFrom: $(runtime.cores)
    - position: 3
      valueFrom: $(runtime.outdir)/MarkedSorted.bam
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
outputs:
    sorted_bam:
        type: File
        outputBinding:
            glob: "MarkedSorted.bam"
        secondaryFiles: [.bai]
