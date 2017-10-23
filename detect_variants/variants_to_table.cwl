#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "SelectVariants (GATK 3.6)"
baseCommand: ["/usr/bin/java", "-Xmx4g", "-jar", "/opt/GenomeAnalysisTK.jar", "-T", "VariantsToTable"]
requirements:
    - class: ResourceRequirement
      ramMin: 6000
      tmpdirMin: 25000
arguments:
    ["-o", { valueFrom: $(runtime.outdir)/variants.tsv }]
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
    fields:
      type:
        type: array
        items: string
        inputBinding:
            prefix: "-F"
      inputBinding:
        position: 3
    genotype_fields:
      type:
        type: array
        items: string
        inputBinding:
            prefix: "-GF"
      inputBinding:
        position: 4
outputs:
    variants_tsv:
        type: File
        outputBinding:
            glob: "variants.tsv"
