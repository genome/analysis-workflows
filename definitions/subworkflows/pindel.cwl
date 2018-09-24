#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "pindel parallel workflow"
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
        secondaryFiles: ["^.crai"]
    normal_cram:
        type: File
        secondaryFiles: ["^.crai"]
    interval_list:
        type: File
    insert_size:
        type: int
        default: 400
outputs:
    unfiltered_vcf:
        type: File
        outputSource: filter/unfiltered_vcf
        secondaryFiles: [".tbi"]
    filtered_vcf:
        type: File
        outputSource: filter/filtered_vcf
        secondaryFiles: [".tbi"]
steps:
    get_chromosome_list:
        run: ../tools/get_chromosome_list.cwl
        in: 
            interval_list: interval_list
        out:
            [chromosome_list]
    pindel_cat:
        scatter: chromosome
        run: pindel_cat.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            chromosome: get_chromosome_list/chromosome_list
            insert_size: insert_size
        out:
            [per_chromosome_pindel_out]
    cat_all:
        run: ../tools/cat_all.cwl
        in:
            chromosome_pindel_outs: pindel_cat/per_chromosome_pindel_out
        out:
            [all_chromosome_pindel_head]
    somaticfilter:
        run: ../tools/pindel_somatic_filter.cwl
        in:
            reference: reference
            pindel_output_summary: cat_all/all_chromosome_pindel_head
        out: 
            [vcf]
    bgzip:
        run: ../tools/bgzip.cwl
        in: 
            file: somaticfilter/vcf
        out:
            [bgzipped_file]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
    region_filter:
        run: ../tools/select_variants.cwl
        in:
            reference: reference
            vcf: index/indexed_vcf
            interval_list: interval_list
        out:
            [filtered_vcf]
    remove_end_tags:
        run: ../tools/remove_end_tags.cwl
        in:
            vcf: region_filter/filtered_vcf
        out:
            [processed_vcf]
    reindex:
        run: ../tools/index_vcf.cwl
        in:
            vcf: remove_end_tags/processed_vcf
        out:
            [indexed_vcf]
    filter:
        run: fp_filter.cwl
        in:
            reference: reference
            cram: tumor_cram
            vcf: reindex/indexed_vcf
            variant_caller: 
                valueFrom: "pindel"
        out:
            [unfiltered_vcf, filtered_vcf]
