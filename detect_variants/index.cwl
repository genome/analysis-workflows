#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "vcf index"
baseCommand: "/usr/bin/tabix"
arguments: ["-p", "vcf"]
requirements:
    - class: InitialWorkDirRequirement
      listing:
          - $(inputs.vcf)
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

