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
      dockerPull: mgibio/cle
arguments:
    ["-genotypeMergeOptions", "PRIORITIZE",
     "--rod_priority_list", "callers,docm",
     "--setKey","null",
     "-o", { valueFrom: $(runtime.outdir)/merged.vcf.gz }]
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "-R"
            position: 1
    callers_vcf:
        type: File
        inputBinding:
            prefix: "--variant:callers"
            position: 2
        secondaryFiles: [.tbi]
    docm_vcf:
        type: File
        inputBinding:
            prefix: "--variant:docm"
            position: 3
outputs:
    merged_vcf:
        type: File
        outputBinding:
            glob: "merged.vcf.gz"
