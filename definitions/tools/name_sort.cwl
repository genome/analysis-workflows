#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'sort BAM by name'
baseCommand: ["/usr/bin/sambamba", "sort"]
arguments: 
    ["-t", { valueFrom: $(runtime.cores) },
    "-m", "12G",
    "-n",
    "-o", { valueFrom: "$(inputs.bam.nameroot).NameSorted.bam" }]
requirements:
    - class: ResourceRequirement
      ramMin: 13000
      coresMin: 8
    - class: DockerRequirement
      dockerPull: "mgibio/sambamba-cwl:0.6.4"
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
outputs:
    name_sorted_bam:
        type: File
        outputBinding:
            glob: "$(inputs.bam.nameroot).NameSorted.bam"
