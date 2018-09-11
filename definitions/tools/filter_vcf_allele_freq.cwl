#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "AF filter"
baseCommand: ["/usr/bin/perl", "/opt/vep/ensembl-vep/filter_vep"]
requirements:
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/cle"
arguments:
    ["--format", "vcf",
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
                            "AF", "<", inputs.maximum_population_allele_frequency,
                            "or not", "AF"
                        ].join(" ")
                    ]
                }
            position: 2
outputs:
    filtered_vcf:
        type: File
        outputBinding:
            glob: "annotated.af_filtered.vcf"
