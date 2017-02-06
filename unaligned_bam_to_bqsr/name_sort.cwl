#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'sort BAM by name'
baseCommand: ["/usr/local/bin/sambamba", "sort"]
arguments: 
    ["-t", "8",
    "-m", "12G",
    "-n",
    "-o", { valueFrom: $(runtime.outdir)/NameSorted.bam }]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/sambamba-0.6.4:1"
    - class: ResourceRequirement
      ramMin: 12000
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
outputs:
    name_sorted_bam:
        type: File
        outputBinding:
            glob: "NameSorted.bam"
