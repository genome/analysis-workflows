#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "SelectVariants (GATK 4.1.8.1)"
baseCommand: ["/gatk/gatk", "--java-options", "-Xmx4g", "SelectVariants"]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.8.1"
arguments:
    ["-O", { valueFrom: $(runtime.outdir)/$(inputs.output_vcf_basename).vcf.gz }]
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
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
            prefix: "--exclude-filtered"
            position: 4
    output_vcf_basename:
        type: string?
        default: select_variants
    samples_to_include:
        type: string[]?
        inputBinding:
            prefix: "--sample-name"
            position: 5
        doc: 'include genotypes from this sample'
    select_type:
        type:
            - "null"
            - type: enum
              symbols: ["INDEL", "SNP", "MIXED", "MNP", "SYMBOLIC", "NO_VARIATION"]
        inputBinding:
            prefix: "-select-type"
            position: 6
        doc: 'select only a certain type of variants' 
outputs:
    filtered_vcf:
        type: File
        secondaryFiles: [.tbi]
        outputBinding:
            glob: $(inputs.output_vcf_basename).vcf.gz
