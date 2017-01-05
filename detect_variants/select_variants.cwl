#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "SelectVariants (GATK 3.6)"
baseCommand: ["/usr/local/bin/jdk1.8.0_45/bin/java", "-jar", "/usr/local/bin/GATK3.6/GenomeAnalysisTK.jar", "-T", "SelectVariants"]
requirements:
    - class: DockerRequirement
      dockerPull: "dbmi/gatk-docker:v1" #GATK 3.6 at a specific container revision
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/output.vcf.gz }]
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
        secondaryFiles: .tbi
    interval_list:
        type: File
        inputBinding:
            prefix: "-L"
            position: 3
outputs:
    filtered_vcf:
        type: File
        outputBinding:
            glob: "output.vcf.gz"
