#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "SelectVariants (GATK 4.1.8.1)"
baseCommand: ["/gatk/gatk", "--java-options", "-Xmx4g", "VariantsToTable"]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
    - class: DockerRequirement
      dockerPull: "broadinstitute/gatk:4.1.8.1"
arguments:
    ["-O", { valueFrom: $(runtime.outdir)/variants.tsv }]
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
    fields:
      type:
        type: array
        items: string
        inputBinding:
            prefix: "-F"
      default: ['CHROM','POS','ID','REF','ALT','set']
      inputBinding:
        position: 3
    genotype_fields:
      type:
        type: array
        items: string
        inputBinding:
            prefix: "-GF"
      default: ['GT','AD','DP','AF']
      inputBinding:
        position: 4
outputs:
    variants_tsv:
        type: File
        outputBinding:
            glob: "variants.tsv"
