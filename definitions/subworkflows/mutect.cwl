#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "mutect parallel workflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    reference:
        type: string
    tumor_cram:
        type: File
        secondaryFiles: [.crai]
    normal_cram:
        type: File?
        secondaryFiles: [.crai]
    interval_list:
        type: File
    scatter_count:
        type: int
    dbsnp_vcf:
        type: File?
        secondaryFiles: [.tbi]
    cosmic_vcf:
        type: File?
        secondaryFiles: [.tbi]
    artifact_detection_mode:
        type: boolean
    panel_of_normals_vcf:
        type: File?
        secondaryFiles: [.tbi]
    max_alt_allele_in_normal_fraction:
        type: float?
    max_alt_alleles_in_normal_count:
        type: int?
outputs:
    unfiltered_vcf:
        type: File
        outputSource: filter/unfiltered_vcf
        secondaryFiles: [.tbi]
    filtered_vcf:
        type: File
        outputSource: filter/filtered_vcf
        secondaryFiles: [.tbi]
steps:
    split_interval_list:
        run: ../tools/split_interval_list.cwl
        in:
            interval_list: interval_list
            scatter_count: scatter_count
        out: [split_interval_lists]
    mutect:
        scatter: interval_list
        run: ../tools/mutect.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            interval_list: split_interval_list/split_interval_lists
            dbsnp_vcf: dbsnp_vcf
            cosmic_vcf: cosmic_vcf
            artifact_detection_mode: artifact_detection_mode
            panel_of_normals_vcf: panel_of_normals_vcf
            max_alt_allele_in_normal_fraction: max_alt_allele_in_normal_fraction
            max_alt_alleles_in_normal_count: max_alt_alleles_in_normal_count
        out:
            [vcf]
    merge:
        run: ../tools/merge_vcf.cwl
        in:
            vcfs: mutect/vcf
        out:
            [merged_vcf]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: merge/merged_vcf
        out:
            [indexed_vcf]
    filter:
        run: fp_filter.cwl
        in:
            reference: reference
            cram: tumor_cram
            vcf: index/indexed_vcf
            variant_caller: 
                valueFrom: "mutect"
        out:
            [unfiltered_vcf, filtered_vcf]


