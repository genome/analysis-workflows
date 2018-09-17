#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Varscan Workflow"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    reference:
        type: string
    cram:
        type: File
        secondaryFiles: [^.crai]
    interval_list:
        type: File
    strand_filter:
        type: int?
        default: 0
    min_coverage:
        type: int?
        default: 8
    min_var_freq:
        type: float?
        default: 0.1
    min_reads:
        type: int?
        default: 2
    p_value:
        type: float?
        default: 0.99
    max_normal_freq:
        type: float?
    sample_name:
        type: string
outputs:
    unfiltered_vcf:
        type: File
        outputSource: filter/unfiltered_vcf
    filtered_vcf:
        type: File
        outputSource: filter/filtered_vcf
        secondaryFiles: [.tbi]
steps:
    intervals_to_bed:
        run: ../tools/intervals_to_bed.cwl
        in:
            interval_list: interval_list
        out:
            [interval_bed]
    varscan:
        run: ../tools/varscan_germline.cwl
        in:
            reference: reference
            cram: cram
            roi_bed: intervals_to_bed/interval_bed
            strand_filter: strand_filter
            min_coverage: min_coverage
            min_var_freq: min_var_freq
            min_reads: min_reads
            p_value: p_value
            sample_name: sample_name
        out:
            [variants]
    bgzip_and_index:
        run: bgzip_and_index.cwl
        in:
            vcf: varscan/variants
        out:
            [indexed_vcf]
    filter:
        run: fp_filter.cwl
        in:
            reference: reference
            cram: cram
            vcf: bgzip_and_index/indexed_vcf
            variant_caller:
                valueFrom: "varscan"
            sample_name: sample_name
            min_var_freq: min_var_freq
        out:
            [unfiltered_vcf, filtered_vcf]
