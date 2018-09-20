#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools index"
arguments: [
    "cp", $(inputs.bam.path), "$(runtime.outdir)/$(inputs.bam.basename)",
    { valueFrom: " && ", shellQuote: false },
    "/opt/samtools/bin/samtools", "index", $(inputs.bam.path), "$(runtime.outdir)/$(inputs.bam.basename).bai"
]
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
inputs:
    bam:
        type: File
outputs:
    indexed_bam:
        type: File
        secondaryFiles: [.bai]
        outputBinding:
            glob: $(inputs.bam.basename)
