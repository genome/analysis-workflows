#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "mark duplicates and sort"

baseCommand: ["/bin/bash", "/usr/bin/markduplicates_helper.sh"]
requirements:
    - class: ResourceRequirement
      coresMin: 8
      ramMin: 18000
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
