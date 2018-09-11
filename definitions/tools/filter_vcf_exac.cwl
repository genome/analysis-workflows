#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "ExAC filter"
baseCommand: ["/usr/bin/perl", "/usr/bin/vcf_check.pl"]
requirements:
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: mgibio/cle
arguments:
    [{ valueFrom: $(inputs.vcf.path) },
    { valueFrom: $(runtime.outdir)/annotated.exac_filtered.vcf },
    "/usr/bin/perl", "/opt/vep/ensembl-vep/filter_vep",
    "--format", "vcf",
    "-o", { valueFrom: $(runtime.outdir)/annotated.exac_filtered.vcf }]
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
                            "(",
                            "ExAC_AF",  "<", inputs.maximum_population_allele_frequency,
                            "and",
                            "ExAC_Adj_AF", "<", inputs.maximum_population_allele_frequency,
                            "and",
                            "ExAC_AFR_AF", "<", inputs.maximum_population_allele_frequency,
                            "and",
                            "ExAC_SAS_AF", "<", inputs.maximum_population_allele_frequency,
                            "and",
                            "ExAC_EAS_AF", "<", inputs.maximum_population_allele_frequency,
                            "and",
                            "ExAC_NFE_AF", "<", inputs.maximum_population_allele_frequency,
                            ")",
                            "or not", "ExAC_AF"
                        ].join(" ")
                    ]
                }
            position: 2
outputs:
    filtered_vcf:
        type: File
        outputBinding:
            glob: "annotated.exac_filtered.vcf"
