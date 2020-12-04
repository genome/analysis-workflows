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
    blocklist_bedpe:
        type: File?
outputs:
    bcftools_merged_sv_vcf:
        type: File
        outputSource: filter_blocklist_bcftools/filtered_sv_vcf
    bcftools_merged_annotated_tsv:
        type: File
        outputSource: bcftools_annotate_variants/sv_variants_tsv
    bcftools_merged_filtered_annotated_tsv:
       type: File
       outputSource: bcftools_annotsv_filter/filtered_tsv
    survivor_merged_sv_vcf:
        type: File
        outputSource: filter_blocklist_survivor/filtered_sv_vcf
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
    filter_blocklist_survivor:
        run: ../tools/filter_sv_vcf_blocklist_bedpe.cwl
        in:
            input_vcf: survivor_merge_sv_vcfs/merged_vcf
            blocklist_bedpe: blocklist_bedpe
            output_vcf_basename:
                default: "SURVIVOR-sv-merged"
        out:
            [filtered_sv_vcf]
    survivor_annotate_variants:
        run: ../tools/annotsv.cwl
        in:
            genome_build: genome_build
            input_vcf: filter_blocklist_survivor/filtered_sv_vcf
            output_tsv_name:
                default: "SURVIVOR-merged-AnnotSV.tsv"
            snps_vcf:
                source: [snps_vcf]
                valueFrom: ${ return [ self ]; }
        out:
            [sv_variants_tsv]
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
    filter_blocklist_bcftools:
        run: ../tools/filter_sv_vcf_blocklist_bedpe.cwl
        in:
            input_vcf: bcftools_merge_sv_vcfs/merged_sv_vcf
            blocklist_bedpe: blocklist_bedpe
            output_vcf_basename:
                default: "bcftools-sv-merged"
        out:
            [filtered_sv_vcf]
    bcftools_annotate_variants:
        run: ../tools/annotsv.cwl
        in:
            genome_build: genome_build
            input_vcf: filter_blocklist_bcftools/filtered_sv_vcf
            output_tsv_name:
                default: "bcftools-merged-AnnotSV.tsv"
            snps_vcf:
                source: [snps_vcf]
                valueFrom: ${ return [ self ]; }
        out:
            [sv_variants_tsv]
    bcftools_annotsv_filter:
        run: ../tools/annotsv_filter.cwl
        in:
            annotsv_tsv: bcftools_annotate_variants/sv_variants_tsv
            filtering_frequency:
                default: "0.05"
        out:
            [filtered_tsv]
