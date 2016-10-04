#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "CombineVariants (GATK 3.6)"
baseCommand: ["/usr/local/bin/jdk1.8.0_45/bin/java", "-jar", "/usr/local/bin/GATK3.6/GenomeAnalysisTK.jar", "-T", "MuTect2"]
requirements:
    - class: DockerRequirement
      dockerPull: "dbmi/gatk-docker:latest@sha256:7e7d7911d51b0109cd1db1b5ac824287c7af013da29fe28022289ff62940f0be" #GATK 3.6 at a specific container revision
arguments:
    ["-genotypeMergeOptions", "PRIORITIZE",
     "--rod_priority_list", "mutect,varscan,strelka",
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
        secondaryFiles: .tbi
    varscan_vcf:
        type: File
        inputBinding:
            prefix: "--variant:varscan"
            position: 3
        secondaryFiles: .tbi
    strelka_vcf:
        type: File
        inputBinding:
            prefix: "--variant:strelka"
            position: 4
        secondaryFiles: .tbi
outputs:
    combined_vcf:
        type: File
        outputBinding:
            glob: "combined.vcf"

