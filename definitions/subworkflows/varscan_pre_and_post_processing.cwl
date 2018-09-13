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
    tumor_cram:
        type: File
        secondaryFiles: [^.crai]
    normal_cram:
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
    p_value:
        type: float?
        default: 0.99
    max_normal_freq:
        type: float?
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
    intervals_to_bed:
        run: ../tools/intervals_to_bed.cwl
        in:
            interval_list: interval_list
        out:
            [interval_bed]
    varscan:
        run: varscan.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            roi_bed: intervals_to_bed/interval_bed
            strand_filter: strand_filter
            min_coverage: min_coverage
            min_var_freq: min_var_freq
            p_value: p_value
            max_normal_freq: max_normal_freq
        out:
            [somatic_snvs, somatic_indels, somatic_hc_snvs, somatic_hc_indels]
    bgzip_and_index_snvs:
        run: bgzip_and_index.cwl
        in:
            vcf: varscan/somatic_snvs
        out:
            [indexed_vcf]
    bgzip_and_index_hc_snvs:
        run: bgzip_and_index.cwl
        in:
            vcf: varscan/somatic_hc_snvs
        out:
            [indexed_vcf]
    bgzip_and_index_indels:
        run: bgzip_and_index.cwl
        in:
            vcf: varscan/somatic_indels
        out:
            [indexed_vcf]
    bgzip_and_index_hc_indels:
        run: bgzip_and_index.cwl
        in:
            vcf: varscan/somatic_hc_indels
        out:
            [indexed_vcf]
    merge_snvs:
        run: ../tools/set_filter_status.cwl
        in:
            vcf: bgzip_and_index_snvs/indexed_vcf
            filtered_vcf: bgzip_and_index_hc_snvs/indexed_vcf
            reference: reference
        out:
            [merged_vcf]
    index_snvs:
        run: ../tools/index_vcf.cwl
        in:
            vcf: merge_snvs/merged_vcf
        out:
            [indexed_vcf]
    merge_indels:
        run: ../tools/set_filter_status.cwl
        in:
            vcf: bgzip_and_index_indels/indexed_vcf
            filtered_vcf: bgzip_and_index_hc_indels/indexed_vcf
            reference: reference
        out:
            [merged_vcf]
    index_indels:
        run: ../tools/index_vcf.cwl
        in:
            vcf: merge_indels/merged_vcf
        out:
            [indexed_vcf]
    merge:
        run: ../tools/merge_vcf.cwl
        in:
            vcfs: [index_snvs/indexed_vcf, index_indels/indexed_vcf]
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
            min_var_freq: min_var_freq
            variant_caller: 
                valueFrom: "varscan"
        out:
            [unfiltered_vcf, filtered_vcf]
