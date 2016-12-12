#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'tag BAM'
baseCommand: "/usr/local/bin/samblaster"
arguments:
    ["-a",
    "--addMateTags",
    "-o", { valueFrom: $(runtime.outdir)/MergedTagged.sam }]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/samblaster-0.1.23:1"
inputs:
    sam:
        type: File
        inputBinding:
            prefix: "-i"
            position: 1
outputs:
    tagged_sam:
        type: File
        outputBinding:
            glob: "MergedTagged.sam"
