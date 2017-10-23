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
    tumor_cram:
        type: File
        secondaryFiles: [.crai]
    normal_cram:
        type: File
        secondaryFiles: [.crai]
    reference:
        type: string
    interval_list:
        type: File
    exome_mode:
        type: boolean
    cpu_reserved:
        type: int?
        default: 8
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
        run: strelka.cwl
        in:
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            reference: reference
            exome_mode: exome_mode
            cpu_reserved: cpu_reserved
        out:
            [indels, snvs]
    process:
        scatter: vcf
        run: process_vcf.cwl
        in:
            vcf: [strelka/snvs, strelka/indels]
        out:
            [processed_vcf]
    merge:
        run: ../detect_variants/merge.cwl
        in:
            vcfs: process/processed_vcf
        out:
            [merged_vcf]
    index_full:
        run: ../detect_variants/index.cwl
        in:
            vcf: merge/merged_vcf
        out:
            [indexed_vcf]
    region_filter:
        run: ../detect_variants/select_variants.cwl
        in:
            reference: reference
            vcf: index_full/indexed_vcf
            interval_list: interval_list
        out:
            [filtered_vcf]
    filter:
        run: ../fp_filter/workflow.cwl
        in:
            reference: reference
            cram: tumor_cram
            vcf: region_filter/filtered_vcf
            variant_caller: 
                valueFrom: "strelka"
        out:
            [unfiltered_vcf, filtered_vcf]



