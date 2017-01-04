#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "verify BAM ID"
baseCommand: "/usr/local/bin/verifyBamID"
arguments:
    ["--out", { valueFrom: $(runtime.outdir)/VerifyBamId }]
requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/genome/lims-verifybamid:1"
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
            glob: "VerifyBamId.selfSM"
    verify_bam_id_depth:
        type: File
        outputBinding:
            glob: "VerifyBamId.depthSM"

