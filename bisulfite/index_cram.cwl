#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "samtools index"
arguments: [
    "cp", $(inputs.cram.path), "$(runtime.outdir)/$(inputs.cram.basename)",
    { valueFrom: " && ", shellQuote: false },
    "/opt/samtools/bin/samtools", "index", $(inputs.cram.path), "$(runtime.outdir)/$(inputs.cram.basename).bai"
]
requirements:
    - class: ShellCommandRequirement
inputs:
    cram:
        type: File
outputs:
    indexed_cram:
        type: File
        secondaryFiles: [.crai]
        outputBinding:
            glob: $(inputs.cram.basename)
