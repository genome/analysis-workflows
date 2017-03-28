#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "CombineVariants (GATK 3.6)"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "CombineVariants"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000

arguments:
    ["-genotypeMergeOptions", "PRIORITIZE",
     "--rod_priority_list", "filtered,docm",
     "-o", { valueFrom: $(runtime.outdir)/combine_docm.vcf }]
inputs:
    reference:
        type: File
        inputBinding:
            prefix: "-R"
            position: 1
        secondaryFiles: [.fai, ^.dict]
    filtered_vcf:
        type: File
        inputBinding:
            prefix: "--variant:filtered"
            position: 2
        secondaryFiles: [.tbi]
    docm_vcf:
        type: File
        inputBinding:
            prefix: "--variant:docm"
            position: 3
        secondaryFiles: [.tbi]
outputs:
    combine_docm_vcf:
        type: File
        outputBinding:
            glob: "combine_docm.vcf"

