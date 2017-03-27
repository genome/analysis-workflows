#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "HaplotypeCaller (GATK 3.6)"
baseCommand: ["/usr/local/bin/jdk1.8.0_45/bin/java", "-jar", "/usr/local/bin/GATK3.6/GenomeAnalysisTK.jar", "-T", "HaplotypeCaller"]
requirements:
    - class: DockerRequirement
      dockerPull: "dbmi/gatk-docker:v1" #GATK 3.6 at a specific container revision
arguments:
    ["-gt_mode", "GENOTYPE_GIVEN_ALLELES",
    "-o", { valueFrom: $(runtime.outdir)/docm_out.vcf }]
inputs:
    reference:
        type: File
        inputBinding:
            prefix: "-R"
            position: 1
        secondaryFiles: [.fai, ^.dict]
    tumor_bam:
        type: File
        inputBinding:
            prefix: "-I"
            position: 2
        secondaryFiles: [^.bai]
    normal_bam:
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
    interval_file:
        type: File
        inputBinding:
            prefix: "-L"
            position: 5
        secondaryFiles: [.tbi]
outputs:
    docm_out:
        type: File
        outputBinding:
            glob: "docm_out.vcf"
