#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Ensembl Variant Effect Predictor"
baseCommand: ["/usr/bin/perl", "/usr/bin/variant_effect_predictor.pl"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 32000
      tmpdirMin: 25000
arguments:
    ["--cache",
    "--offline",
    "--format", "vcf",
    "--vcf",
    "--plugin", "Downstream",
    "--plugin", "Wildtype",
    "--symbol",
    "--term", "SO",
    "--flag_pick",
    "--maf_exac",
    "-o", { valueFrom: $(runtime.outdir)/annotated.vcf },
    "--dir",
    { valueFrom: "`", shellQuote: false },
    { valueFrom: "cat", shellQuote: false },
    { valueFrom: $(inputs.cache_dir), shellQuote: false },
    { valueFrom: "`", shellQuote: false }]
inputs:
    vcf:
        type: File
        inputBinding:
            prefix: "-i"
            position: 1
    cache_dir:
        type: File
    synonyms_file:
        type: File?
        inputBinding:
            prefix: "--synonyms"
            position: 2
outputs:
    annotated_vcf:
        type: File
        outputBinding:
            glob: "annotated.vcf"
