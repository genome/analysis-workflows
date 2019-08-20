#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

baseCommand: ["/opt/bcftools/bin/bcftools", "merge"]

requirements:
    - class: ResourceRequirement
      ramMin: 4000
    - class: DockerRequirement
      dockerPull: "mgibio/bcftools-cwl:1.9"

inputs:
    force_merge:
        type: boolean?
        default: true
        inputBinding:
            position: 1
            prefix: "--force-samples"
        doc: "resolve duplicate sample names"
    merge_method:
        type:
            type: enum
            symbols: ["none", "snps", "indels", "both", "all", "id"]
        default: "none"
        inputBinding:
            position: 2
            prefix: "--merge"
        doc: "method used to merge allow multiallelic indels/snps"
    missing_ref:
        type: boolean?
        default: false
        inputBinding:
            position: 3
            prefix: "--missing-to-ref"
        doc: "assume genotypes at missing sites are 0/0"
    output_type:
        type:
            type: enum
            symbols: ["b", "u", "z", "v"]
        default: "z"
        inputBinding:
            position: 4
            prefix: "--output-type"
        doc: "output file format"
    output_vcf_name:
        type: string?
        default: "bcftools_merged.vcf.gz"
        inputBinding:
            position: 5
            prefix: "--output"
        doc: "output vcf file name"
    vcfs:
        type: File[]
        inputBinding:
            position: 6
        doc: "input bgzipped tabix indexed vcfs to merge"

outputs:
    merged_sv_vcf:
        type: File
        outputBinding:
            glob: $(inputs.output_vcf_name)

