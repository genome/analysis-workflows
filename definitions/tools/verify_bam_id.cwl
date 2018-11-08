#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "verify BAM ID"
baseCommand: "/usr/local/bin/verifyBamID"
arguments:
    ["--out", { valueFrom: $(runtime.outdir)/$(inputs.bam.nameroot).VerifyBamId }]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/verify_bam_id-cwl:1.1.3"
inputs:
    vcf:
        type: File
        inputBinding:
            prefix: "--vcf"
    bam:
        type: File
        inputBinding:
            prefix: "--bam"
outputs:
    verify_bam_id_metrics:
        type: File
        outputBinding:
            glob: "$(inputs.bam.nameroot).VerifyBamId.selfSM"
    verify_bam_id_depth:
        type: File
        outputBinding:
            glob: "$(inputs.bam.nameroot).VerifyBamId.depthSM"

