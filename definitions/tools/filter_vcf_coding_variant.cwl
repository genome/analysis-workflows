#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Coding Variant filter"
baseCommand: ["/usr/bin/perl", "/usr/bin/vcf_check.pl"]
requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/vep_helper-cwl:1.0.0"
arguments:
    [{ valueFrom: $(inputs.vcf.path) },
    { valueFrom: $(runtime.outdir)/annotated.coding_variant_filtered.vcf },
    "/usr/bin/perl", "/opt/vep/src/ensembl-vep/filter_vep",
    "--format", "vcf",
    "-o", { valueFrom: $(runtime.outdir)/annotated.coding_variant_filtered.vcf },
    "--ontology",
    "--filter", "Consequence is coding_sequence_variant"]
inputs:
    vcf:
        type: File
        inputBinding:
            prefix: "-i"
            position: 1
outputs:
    filtered_vcf:
        type: File
        outputBinding:
            glob: "annotated.coding_variant_filtered.vcf"
