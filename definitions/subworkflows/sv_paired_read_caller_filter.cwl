#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Filter single sample sv vcf from paired read callers(Manta/Smoove)"

inputs:
    abundance_percentage:
        type: double?
    bam:
        type: File
    deletion_depth:
        type: double?
    duplication_depth:
        type: double?
    output_vcf_name:
        type: string?
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    snps_vcf:
        type: File?
    sv_paired_count:
        type: int?
    sv_split_count:
        type: int?
    sv_vcf:
        type: File
    vcf_source:
        type:
          - type: enum
            symbols: ["manta", "smoove"]
outputs:
    filtered_vcf:
        type: File
        outputSource: filtered_vcf_index/indexed_vcf
steps:
    read_support_filter:
        run: ../tools/filter_sv_vcf_read_support.cwl
        in:
            abundance_percentage: abundance_percentage
            input_vcf: sv_vcf
            paired_count: sv_paired_count
            split_count: sv_split_count
            vcf_source: vcf_source
        out:
            [filtered_sv_vcf]
    duphold_annotate:
        run: ../tools/duphold.cwl
        in:
            bam: bam
            reference: reference
            snps_vcf: snps_vcf
            sv_vcf: read_support_filter/filtered_sv_vcf
        out:
            [annotated_sv_vcf]
    depth_filter:
        run: ../tools/filter_sv_vcf_depth.cwl
        in:
            input_vcf: duphold_annotate/annotated_sv_vcf
            deletion_depth: deletion_depth
            duplication_depth: duplication_depth
            output_vcf_name: output_vcf_name
            vcf_source:
                default: "duphold"
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
