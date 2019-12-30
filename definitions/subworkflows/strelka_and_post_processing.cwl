#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "strelka workflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    tumor_bam:
        type: File
        secondaryFiles: [.bai,^.bai]
    normal_bam:
        type: File
        secondaryFiles: [.bai,^.bai]
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    interval_list:
        type: File
    exome_mode:
        type: boolean
    cpu_reserved:
        type: int?
        default: 8
    normal_sample_name:
        type: string
    tumor_sample_name:
        type: string
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
    strelka:
        run: ../tools/strelka.cwl
        in:
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            reference: reference
            exome_mode: exome_mode
            cpu_reserved: cpu_reserved
        out:
            [indels, snvs]
    process:
        scatter: vcf
        run: strelka_process_vcf.cwl
        in:
            vcf: [strelka/snvs, strelka/indels]
        out:
            [processed_vcf]
    merge:
        run: ../tools/merge_vcf.cwl
        in:
            vcfs: process/processed_vcf
        out:
            [merged_vcf]
    rename_tumor_sample:
        run: ../tools/replace_vcf_sample_name.cwl
        in:
            input_vcf: merge/merged_vcf
            sample_to_replace:
                default: 'TUMOR'
            new_sample_name: tumor_sample_name
        out:
            [renamed_vcf]
    rename_normal_sample:
        run: ../tools/replace_vcf_sample_name.cwl
        in:
            input_vcf: rename_tumor_sample/renamed_vcf
            sample_to_replace:
                default: 'NORMAL'
            new_sample_name: normal_sample_name
        out:
            [renamed_vcf]
    index_full:
        run: ../tools/index_vcf.cwl
        in:
            vcf: rename_normal_sample/renamed_vcf
        out:
            [indexed_vcf]
    region_filter:
        run: ../tools/select_variants.cwl
        in:
            reference: reference
            vcf: index_full/indexed_vcf
            interval_list: interval_list
        out:
            [filtered_vcf]
    filter:
        run: fp_filter.cwl
        in:
            reference: reference
            bam: tumor_bam
            vcf: region_filter/filtered_vcf
            variant_caller: 
                valueFrom: "strelka"
            sample_name: tumor_sample_name
        out:
            [unfiltered_vcf, filtered_vcf]
