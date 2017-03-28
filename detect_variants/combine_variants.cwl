#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "CombineVariants (GATK 3.6)"
<<<<<<< HEAD
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "CombineVariants"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
      tmpdirMin: 25000
arguments:
    ["-genotypeMergeOptions", "PRIORITIZE",
     "--rod_priority_list", "mutect,varscan,strelka,pindel",
     "-o", { valueFrom: $(runtime.outdir)/combined.vcf }]
inputs:
    reference:
        type: File
        inputBinding:
            prefix: "-R"
            position: 1
        secondaryFiles: [".fai", "^.dict"]
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
outputs:
    combined_vcf:
        type: File
        outputBinding:
            glob: "combined.vcf"

