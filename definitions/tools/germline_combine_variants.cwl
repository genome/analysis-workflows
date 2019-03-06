#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "CombineVariants (GATK 3.6)"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "CombineVariants"]
requirements:
    - class: ResourceRequirement
      ramMin: 9000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "mgibio/gatk-cwl:3.6.0"
arguments:
    ["-genotypeMergeOptions", "PRIORITIZE",
     "--rod_priority_list", "varscan,docm",
     "-o", { valueFrom: $(runtime.outdir)/combined.vcf.gz }]
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "-R"
            position: 1
    varscan_vcf:
        type: File
        inputBinding:
            prefix: "--variant:varscan"
            position: 3
        secondaryFiles: [.tbi]
    docm_vcf:
        type: File
        inputBinding:
            prefix: "--variant:docm"
            position: 6
        secondaryFiles: [.tbi]
outputs:
    combined_vcf:
        type: File
        outputBinding:
            glob: "combined.vcf.gz"

