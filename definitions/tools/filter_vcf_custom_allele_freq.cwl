#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Custom allele frequency filter"
baseCommand: ["/usr/bin/perl", "/usr/bin/vcf_check.pl"]
requirements:
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/vep_helper-cwl:1.1.0"
    - class: ResourceRequirement
      ramMin: 4000
    - class: StepInputExpressionRequirement

arguments:
    [{ valueFrom: $(inputs.vcf.path) },
    { valueFrom: $(runtime.outdir)/annotated.af_filtered.vcf },
    "/usr/bin/perl", "/opt/vep/src/ensembl-vep/filter_vep",
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
                            inputs.field_name, "<", inputs.maximum_population_allele_frequency,
                            "or not", inputs.field_name
                        ].join(" ")
                    ]
                }
            position: 2
    field_name:
        type: string
outputs:
    filtered_vcf:
        type: File
        outputBinding:
            glob: "annotated.af_filtered.vcf"
