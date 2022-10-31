#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Normalize variants"
baseCommand: ["/gatk/gatk", "--java-options", "-Xmx8g", "LeftAlignAndTrimVariants"]
requirements:
    - class: ResourceRequirement
      ramMin: 9000
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.8.1"
    - class: InitialWorkDirRequirement
      listing:
      - $(inputs.ref_fai)
      - $(inputs.ref_dict)
      - $(inputs.reference)

arguments:
    ["-O", { valueFrom: $(runtime.outdir)/normalized.vcf.gz }]
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
        inputBinding:
            prefix: "-R"
            position: 1
    ref_fai:
        type: File
    ref_dict:
        type: File
    vcf:
        type: File
        inputBinding:
            prefix: "-V"
            position: 2
        secondaryFiles: [".tbi"]
outputs:
    normalized_vcf:
        type: File
        secondaryFiles: [".tbi"]
        outputBinding:
            glob: "normalized.vcf.gz"
