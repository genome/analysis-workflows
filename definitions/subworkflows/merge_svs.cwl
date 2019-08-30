#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Merge, annotate, and generate a TSV for SVs"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement

inputs:
    cohort_name:
        type: string?
    estimate_sv_distance:
        type: boolean
    genome_build:
        type:
            type: enum
            symbols: ['GRCh37', 'GRCh38', 'GRCm37', 'GRCm38']
    max_distance_to_merge:
        type: int
    minimum_sv_calls:
        type: int
    minimum_sv_size:
        type: int
    same_strand:
        type: boolean
    same_type:
        type: boolean
    snps_vcf:
        type: File?
    sv_vcfs:
        type: File[]
outputs:
    merged_sv_vcf:
        type: File
        outputSource: bgzip_merged_sv_vcf/bgzipped_file
    merged_annotated_tsv:
        type: File
        outputSource: annotate_variants/sv_variants_tsv
steps:
    merge_sv_vcfs:
        run: ../tools/survivor.cwl
        in: 
            vcfs: sv_vcfs
            max_distance_to_merge: max_distance_to_merge
            minimum_sv_calls: minimum_sv_calls
            same_type: same_type
            same_strand: same_strand
            estimate_sv_distance: estimate_sv_distance
            minimum_sv_size: minimum_sv_size
            cohort_name: cohort_name
        out:
            [merged_vcf]
    annotate_variants:
        run: ../tools/annotsv.cwl
        in:
            genome_build: genome_build
            input_vcf: merge_sv_vcfs/merged_vcf
            snps_vcf:
                source: [snps_vcf]
                valueFrom: ${ return [ self ]; }
        out:
            [sv_variants_tsv]
    bgzip_merged_sv_vcf:
        run: ../tools/bgzip.cwl
        in:
            file: merge_sv_vcfs/merged_vcf
        out:
            [bgzipped_file]
