#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "vcf index"
requirements:
    - class: DockerRequirement
      dockerPull: mgibio/cle
arguments: [
    "cp", $(inputs.vcf.path), "$(runtime.outdir)/$(inputs.vcf.basename)",
    { valueFrom: " && ", shellQuote: false },
    "/usr/bin/tabix", "-p", "vcf"
]
requirements:
    - class: ShellCommandRequirement
inputs:
    vcf:
        type: File
        inputBinding:
            valueFrom:
                $(self.basename)
            position: 1
outputs:
    indexed_vcf:
        type: File
        secondaryFiles: [.tbi]
        outputBinding:
            glob: $(inputs.vcf.basename)

