#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "gnomADe_AF filter"
baseCommand: ["/usr/bin/perl", "/usr/bin/vcf_check.pl"]
requirements:
    - class: InlineJavascriptRequirement
arguments:
    [{ valueFrom: $(inputs.vcf.path) },
    { valueFrom: $(runtime.outdir)/annotated.af_filtered.vcf },
    "/usr/bin/perl", "/opt/vep/ensembl-vep/filter_vep",
    "--format", "vcf",
    "-o", { valueFrom: $(runtime.outdir)/annotated.af_filtered.vcf }]
inputs:
    vcf:
        type: File
        inputBinding:
            prefix: "-i"
            position: 1
    maximum_population_allele_frequency:
        type: float
        inputBinding:
            valueFrom: |
                ${
                    return [
                        "--filter",
                        [
                            "gnomADe_AF", "<", inputs.maximum_population_allele_frequency,
                            "or not", "gnomADe_AF"
                        ].join(" ")
                    ]
                }
            position: 2
outputs:
    gnomADe_allele_freq_filtered_vcf:
        type: File
        outputBinding:
            glob: "annotated.af_filtered.vcf"
