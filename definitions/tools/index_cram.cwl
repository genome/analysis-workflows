#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: 'samtools index cram'
arguments: [
    "/opt/samtools/bin/samtools", "index", "$(runtime.outdir)/$(inputs.cram.basename)", "$(runtime.outdir)/$(inputs.cram.basename).crai",
    { valueFrom: " && ", shellQuote: false },
    "ln", "-s", "$(inputs.cram.basename).crai", "$(runtime.outdir)/$(inputs.cram.nameroot).crai"
]
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
    - class: ResourceRequirement
      ramMin: 4000
    - class: InitialWorkDirRequirement
      listing:
          - $(inputs.cram)
inputs:
    cram:
        type: File
outputs:
    indexed_cram:
        type: File
        secondaryFiles: [.crai, ^.crai]
        outputBinding:
            glob: $(inputs.cram.basename)
