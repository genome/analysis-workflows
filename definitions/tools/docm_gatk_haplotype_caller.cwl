#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "HaplotypeCaller (GATK 3.6)"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "HaplotypeCaller"]
requirements:
    - class: ResourceRequirement
      ramMin: 9000
    - class: DockerRequirement
      dockerPull: "mgibio/gatk-cwl:3.6.0"
arguments:
    ["-gt_mode", "GENOTYPE_GIVEN_ALLELES",
    "-o", { valueFrom: $(runtime.outdir)/docm_raw_variants.vcf }]
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "-R"
            position: 1
    normal_bam:
        type: File?
        inputBinding:
            prefix: "-I"
            position: 2
        secondaryFiles: [^.bai]
    bam:
        type: File
        inputBinding:
            prefix: "-I"
            position: 3
        secondaryFiles: [^.bai]
    docm_vcf:
        type: File
        inputBinding:
            prefix: "--alleles"
            position: 4
        secondaryFiles: [.tbi]
    interval_list:
        type: File
        inputBinding:
            prefix: "-L"
            position: 5
outputs:
    docm_raw_variants:
        type: File
        outputBinding:
            glob: "docm_raw_variants.vcf"
