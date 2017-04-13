#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Ensembl Variant Effect Predictor"
baseCommand: ["/usr/bin/perl", "-I", "/opt/lib/perl/VEP/Plugins", "/usr/bin/variant_effect_predictor.pl"]
requirements:
    - class: ShellCommandRequirement
    - class: InlineJavascriptRequirement
    - class: ResourceRequirement
      ramMin: 32000
      tmpdirMin: 25000
arguments:
    ["--format", "vcf",
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
        type: string?
        inputBinding:
            valueFrom: |
                ${
                    if (inputs.cache_dir) {
                        return ["--offline", "--cache", "--maf_exac", "--dir", inputs.cache_dir ]
                    }
                    else {
                        return "--database"
                    }
                }
            position: 4
    synonyms_file:
        type: File?
        inputBinding:
            prefix: "--synonyms"
            position: 2
    coding_only:
        type: boolean
        inputBinding:
            prefix: "--coding_only"
            position: 3
        default: false
outputs:
    annotated_vcf:
        type: File
        outputBinding:
            glob: "annotated.vcf"
    vep_summary:
        type: File
        outputBinding:
            glob: "annotated.vcf_summary.html"
