#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Merge, annotate, and generate a TSV for SVs"
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement

inputs:
    cohort_name:
        type: string?
    estimate_sv_distance:
        type: boolean
    genome_build:
        type: string
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
    bcftools_merged_sv_vcf:
        type: File
        outputSource: bcftools_bgzip_merged_sv_vcf/bgzipped_file
    bcftools_merged_annotated_tsv:
        type: File
        outputSource: bcftools_annotate_variants/sv_variants_tsv
    survivor_merged_sv_vcf:
        type: File
        outputSource: survivor_bgzip_merged_sv_vcf/bgzipped_file
    survivor_merged_annotated_tsv:
        type: File
        outputSource: survivor_annotate_variants/sv_variants_tsv
steps:
    survivor_merge_sv_vcfs:
        run: ../tools/survivor.cwl
        in: 
            vcfs: sv_vcfs
            max_distance_to_merge: max_distance_to_merge
            minimum_sv_calls: minimum_sv_calls
            same_type: same_type
            same_strand: same_strand
            estimate_sv_distance: estimate_sv_distance
            minimum_sv_size: minimum_sv_size
            cohort_name:
                default: "SURVIVOR-sv-merged.vcf"
        out:
            [merged_vcf]
    survivor_annotate_variants:
        run: ../tools/annotsv.cwl
        in:
            genome_build: genome_build
            input_vcf: survivor_merge_sv_vcfs/merged_vcf
            output_tsv_name:
                default: "SURVIVOR-merged-AnnotSV.tsv"
            snps_vcf:
                source: [snps_vcf]
                valueFrom: ${ return [ self ]; }
        out:
            [sv_variants_tsv]
    survivor_bgzip_merged_sv_vcf:
        run: ../tools/bgzip.cwl
        in:
            file: survivor_merge_sv_vcfs/merged_vcf
        out:
            [bgzipped_file]
    bcftools_merge_sv_vcfs:
        run: ../tools/bcftools_merge.cwl
        in:
            merge_method:
                default: "none"
            output_type:
                default: "v"
            output_vcf_name:
                default: "bcftools-sv-merged.vcf"
            vcfs: sv_vcfs
        out:
            [merged_sv_vcf]
    bcftools_annotate_variants:
        run: ../tools/annotsv.cwl
        in:
            genome_build: genome_build
            input_vcf: bcftools_merge_sv_vcfs/merged_sv_vcf
            output_tsv_name:
                default: "bcftools-merged-AnnotSV.tsv"
            snps_vcf:
                source: [snps_vcf]
                valueFrom: ${ return [ self ]; }
        out:
            [sv_variants_tsv]
    bcftools_bgzip_merged_sv_vcf:
        run: ../tools/bgzip.cwl
        in:
            file: bcftools_merge_sv_vcfs/merged_sv_vcf
        out:
            [bgzipped_file]
