#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Normalize variants"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "LeftAlignAndTrimVariants"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/normalized.vcf.gz }]
inputs:
    reference:
        type: File
        inputBinding:
            prefix: "-R"
            position: 1
        secondaryFiles: [".fai", "^.dict"]
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
