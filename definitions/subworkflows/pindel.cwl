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
    tumor_bam:
        type: File
        secondaryFiles: [^.bai]
    normal_bam:
        type: File
        secondaryFiles: [^.bai]
    interval_list:
        type: File
    insert_size:
        type: int
        default: 400
    scatter_count:
        type: int
        default: 50
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
    split_interval_list_to_bed:
        run: ../tools/split_interval_list_to_bed.cwl
        in: 
            interval_list: interval_list
            scatter_count: scatter_count
        out:
            [split_beds]
    pindel_cat:
        scatter: region_file
        run: pindel_cat.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            region_file: split_interval_list_to_bed/split_beds
            insert_size: insert_size
        out:
            [per_region_pindel_out]
    cat_all:
        run: ../tools/cat_all.cwl
        in:
            region_pindel_outs: pindel_cat/per_region_pindel_out
        out:
            [all_region_pindel_head]
    somaticfilter:
        run: ../tools/pindel_somatic_filter.cwl
        in:
            reference: reference
            pindel_output_summary: cat_all/all_region_pindel_head
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
    remove_end_tags:
        run: ../tools/remove_end_tags.cwl
        in:
            vcf: index/indexed_vcf
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
            bam: tumor_bam
            vcf: reindex/indexed_vcf
            variant_caller: 
                valueFrom: "pindel"
        out:
            [unfiltered_vcf, filtered_vcf]
