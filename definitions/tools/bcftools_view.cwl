#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/opt/bcftools/bin/bcftools", "view"]

requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.12"

inputs:
    sample_name:
        type: string?
        inputBinding:
            position: 1
            prefix: "--samples"
        doc: "comma-separated list of samples to include (or exclude with '^' prefix)"
    output_type:
        type:
            type: enum
            symbols: ["b", "u", "z", "v"]
        default: "z"
        inputBinding:
            position: 4
            prefix: "--output-type"
    output_vcf_name:
        type: string?
        default: "bcftools_split.vcf.gz"
        inputBinding:
            position: 5
            prefix: "--output-file"
        doc: "output vcf file name"
    variant_type:
        type: string?
        inputBinding:
             position: 6
             prefix: "--types"
        doc: "select comma-separated list of variant types: snps,indels,mnps,ref,bnd,other"
    in_vcf:
        type: File
        inputBinding:
            position: 7
        doc: "input bgzipped tabix indexed vcf to view"

outputs:
    vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)
