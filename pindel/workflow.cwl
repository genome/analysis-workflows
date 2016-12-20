#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "pindel parallel workflow"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type: File
        secondaryFiles: [".fai"]
    tumor_bam:
        type: File
        secondaryFiles: .bai
    normal_bam:
        type: File
        secondaryFiles: .bai
    interval_list:
        type: File
    scatter_count:
        type: int
        default: 50
    insert_size:
        type: int
        default: 400
outputs:
    merged_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: .tbi
steps:
    split_interval_list:
        run: ../detect_variants/split_interval_list.cwl
        in:
            interval_list: interval_list
            scatter_count: scatter_count
        out:
            [split_interval_lists]
    pindel_cat_grep:
        scatter: interval_list
        run: pindel_cat_grep.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: split_interval_list/split_interval_lists
            insert_size: insert_size
        out:
            [per_interval_pindel_head]
    cat_head:
        run: cat_head.cwl
        in:
            interval_pindel_heads: [pindel_cat_grep/per_interval_pindel_head]
        out:
            [all_interval_pindel_head]
    somaticfilter:
        run: somaticfilter.cwl
        in:
            reference: reference
            pindel_output_summary: cat_head/all_interval_pindel_head
        out: 
            [vcf]
    bgzip:
        run: ../bgzip.cwl
        in: 
            file: somaticfilter/vcf
        out:
            [bgzipped_file]
    index:
        run: ../index.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
