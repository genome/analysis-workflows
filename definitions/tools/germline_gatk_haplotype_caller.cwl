#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "HaplotypeCaller (GATK 3.6)"
baseCommand: ["/usr/bin/java", "-Xmx8g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "HaplotypeCaller"]
requirements:
    - class: ResourceRequirement
      ramMin: 8000
    - class: DockerRequirement
      dockerPull: "mgibio/cle"
arguments:
    ["-gt_mode", "GENOTYPE_GIVEN_ALLELES",
    "-o", { valueFrom: $(runtime.outdir)/docm_out.vcf }]
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "-R"
            position: 1
    cram:
        type: File
        inputBinding:
            prefix: "-I"
            position: 2
        secondaryFiles: [^.crai]
    docm_vcf:
        type: File
        inputBinding:
            prefix: "--alleles"
            position: 3
        secondaryFiles: [.tbi]
    interval_list:
        type: File
        inputBinding:
            prefix: "-L"
            position: 4
outputs:
    docm_out:
        type: File
        outputBinding:
            glob: "docm_out.vcf"
