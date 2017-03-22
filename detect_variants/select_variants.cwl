#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "SelectVariants (GATK 3.6)"
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "SelectVariants"]
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/output.vcf.gz }]
requirements:
    - class: ResourceRequirement
      ramMin: 4000
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
        secondaryFiles: [.tbi]
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
