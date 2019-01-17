#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Normalize variants"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "LeftAlignAndTrimVariants"]
requirements:
    - class: ResourceRequirement
      ramMin: 9000
    - class: DockerRequirement
      dockerPull: "mgibio/gatk-cwl:3.6.0"
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/normalized.vcf.gz }]
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "-R"
            position: 1
    vcf:
        type: File
        inputBinding:
            prefix: "--variant"
            position: 2
        secondaryFiles: [".tbi"]
outputs:
    normalized_vcf:
        type: File
        secondaryFiles: [".tbi"]
        outputBinding:
            glob: "normalized.vcf.gz"
