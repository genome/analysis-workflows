#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "run vt normalize"
baseCommand: ["vt", "normalize"]
requirements:
    - class: DockerRequirement
      dockerPull: quay.io/biocontainers/vt:0.57721--hf74b74d_1
    - class: ResourceRequirement
      ramMin: 4000
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/normalized.vcf.gz }]
inputs:
    vcf:
        type: File
        inputBinding:
            position: 1
        secondaryFiles: [".tbi"]
    reference:
        type:
            - string
            - File
        secondaryFiles: [".fai"]
        inputBinding:
            prefix: "-r"
            position: 2
outputs:
    normalized_vcf:
        type: File
        outputBinding:
            glob: "normalized.vcf.gz"
