#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Filter single sample sv vcf from depth callers(cnvkit/cnvnator)"

inputs:
    deletion_depth:
        type: double?
    duplication_depth:
        type: double?
    min_sv_size:
        type: int?
    output_vcf_name:
        type: string?
    sv_vcf:
        type: File
    vcf_source:
        type:
          - type: enum
            symbols: ["cnvkit", "cnvnator"]
outputs:
    filtered_vcf:
        type: File
        outputSource: filtered_vcf_index/indexed_vcf
steps:
    size_filter:
        run: ../tools/filter_sv_vcf_size.cwl
        in:
            input_vcf: sv_vcf
            size_method:
                default: "min_len"
            sv_size: min_sv_size
        out:
            [filtered_sv_vcf]
    depth_filter:
        run: ../tools/filter_sv_vcf_depth.cwl
        in:
            input_vcf: size_filter/filtered_sv_vcf
            deletion_depth: deletion_depth
            duplication_depth: duplication_depth
            output_vcf_name: output_vcf_name
            vcf_source: vcf_source
        out:
            [filtered_sv_vcf]
    filtered_vcf_bgzip:
        run: ../tools/bgzip.cwl
        in:
            file: depth_filter/filtered_sv_vcf
        out: [bgzipped_file]
    filtered_vcf_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: filtered_vcf_bgzip/bgzipped_file
        out:
            [indexed_vcf]
