#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "varscan v2.4.2 processSomatic"
baseCommand: "processSomatic"
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/varscan:v2.4.2"
    - class: InlineJavascriptRequirement
    - class: InitialWorkDirRequirement
      listing:
        - $(inputs.variants)
inputs:
    variants:
        type: File
        inputBinding:
            valueFrom:
                $(runtime.outdir + "/" + self.basename)
            position: 1
outputs:
    somatic_hc:
        type: File
        outputBinding:
            glob: "*.Somatic.hc.vcf"
    somatic:
        type: File
        outputBinding:
            glob: "*.Somatic.vcf"
    germline_hc:
        type: File
        outputBinding:
            glob: "*.Germline.hc.vcf"
    germline:
        type: File
        outputBinding:
            glob: "*.Germline.vcf"
    loh_hc:
        type: File
        outputBinding:
            glob: "*.LOH.hc.vcf"
    loh:
        type: File
        outputBinding:
            glob: "*.LOH.vcf"
