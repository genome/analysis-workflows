#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "POPFREQ_AF filter"
baseCommand: ["/usr/bin/perl", "/usr/bin/vcf_check.pl"]
requirements:
    - class: InlineJavascriptRequirement
    - class: DockerRequirement
      dockerPull: "mgibio/vep_helper-cwl:1.0.1"
    - class: ResourceRequirement
      ramMin: 4000
arguments:
    [{ valueFrom: $(inputs.vcf.path) },
    { valueFrom: $(runtime.outdir)/annotated.af_filtered.vcf },
    "/usr/bin/perl", "/opt/vep/src/ensembl-vep/filter_vep",
    "--format", "vcf",
    "--vcf_info_field", "tmp",
    "-o", { valueFrom: $(runtime.outdir)/annotated.max_sv_pf_filtered.vcf }]
inputs:
    vcf:
        type: File
        inputBinding:
            prefix: "-i"
            position: 1
    maximum_sv_population_frequency:
        type: float
        inputBinding:
            valueFrom: |
                ${
                    return [
                        "--filter",
                        [
                            "POPFREQ_AF", "<", inputs.maximum_sv_population_frequency,
                            "or not", "POPFREQ_AF"
                        ].join(" ")
                    ]
                }
            position: 2
outputs:
    filtered_vcf:
        type: File
        outputBinding:
            glob: "annotated.max_sv_pf_filtered.vcf"
