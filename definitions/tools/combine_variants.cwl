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
     "--rod_priority_list", "mutect,varscan,strelka,pindel,docm",
     "-o", { valueFrom: $(runtime.outdir)/combined.vcf.gz }]
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "-R"
            position: 1
    mutect_vcf:
        type: File
        inputBinding:
            prefix: "--variant:mutect"
            position: 2
        secondaryFiles: [.tbi]
    varscan_vcf:
        type: File
        inputBinding:
            prefix: "--variant:varscan"
            position: 3
        secondaryFiles: [.tbi]
    strelka_vcf:
        type: File
        inputBinding:
            prefix: "--variant:strelka"
            position: 4
        secondaryFiles: [.tbi]
    pindel_vcf:
        type: File
        inputBinding:
            prefix: "--variant:pindel"
            position: 5
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

