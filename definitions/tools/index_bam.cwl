#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools index"
arguments: [
    "/opt/samtools/bin/samtools", "index", "$(runtime.outdir)/$(inputs.bam.basename)", "$(runtime.outdir)/$(inputs.bam.basename).bai",
    { valueFrom: " && ", shellQuote: false },
    "cp", "$(inputs.bam.basename).bai", "$(runtime.outdir)/$(inputs.bam.nameroot).bai"
]
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
    - class: ResourceRequirement
      ramMin: 4000
    - class: InitialWorkDirRequirement
      listing:
        - ${ var f = inputs.bam; delete f.secondaryFiles; return f }
    - class: InlineJavascriptRequirement
inputs:
    bam:
        type: File
outputs:
    indexed_bam:
        type: File
        secondaryFiles: [.bai, ^.bai]
        outputBinding:
            glob: $(inputs.bam.basename)
