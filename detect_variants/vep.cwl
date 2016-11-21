#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Ensembl Variant Effect Predictor"
baseCommand: ["/usr/bin/perl", "/usr/bin/variant_effect_predictor.pl"]

requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/vep-cwl:v86"
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
    "-o", { valueFrom: $(runtime.outdir)/annotated.vcf }]
inputs:
    vcf:
        type: File
        inputBinding:
            prefix: "-i"
            position: 1
    cache_dir:
        type: Directory
        inputBinding:
            prefix: "--dir"
            position: 2
outputs:
    annotated_vcf:
        type: File
        outputBinding:
            glob: "annotated.vcf"
