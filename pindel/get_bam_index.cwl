#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools index"
baseCommand: ["/usr/local/bin/samtools", "index"]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/samtools-1.3.1-2:2"
arguments:
    - position: 2
      valueFrom: $(runtime.outdir)/$(inputs.bam.basename).bai
inputs:
    bam:
        type: File
        inputBinding:
            position: 1
outputs:
    bam_index:
        type: File
        outputBinding:
            glob: $(inputs.bam.basename).bai
