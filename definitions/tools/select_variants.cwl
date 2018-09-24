#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "SelectVariants (GATK 3.6)"
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "SelectVariants"]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "mgibio/gatk-cwl:3.6.0"
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/$(inputs.output_vcf_basename).vcf.gz }]
inputs:
    reference:
        type: string
        inputBinding:
            prefix: "-R"
            position: 1
    vcf:
        type: File
        inputBinding:
            prefix: "--variant"
            position: 2
        secondaryFiles: [.tbi]
    interval_list:
        type: File?
        inputBinding:
            prefix: "-L"
            position: 3
    exclude_filtered:
        type: boolean?
        inputBinding:
            prefix: "--excludeFiltered"
            position: 4
    output_vcf_basename:
        type: string?
        default: select_variants
outputs:
    filtered_vcf:
        type: File
        secondaryFiles: [.tbi]
        outputBinding:
            glob: $(inputs.output_vcf_basename).vcf.gz
