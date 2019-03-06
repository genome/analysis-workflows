#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "mark duplicates and sort"

baseCommand: ["/bin/bash", "/usr/bin/markduplicates_helper.sh"]
requirements:
    - class: ResourceRequirement
      coresMin: 8
      ramMin: 40000
    - class: DockerRequirement
      dockerPull: "mgibio/mark_duplicates-cwl:1.0.1"
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
