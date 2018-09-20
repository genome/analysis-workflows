#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "vcf index"
arguments: [
    "cp", $(inputs.vcf.path), "$(runtime.outdir)/$(inputs.vcf.basename)",
    { valueFrom: " && ", shellQuote: false },
    "/usr/bin/tabix", "-p", "vcf"
]
requirements:
    - class: ShellCommandRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/samtools-cwl:1.0.0"
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

