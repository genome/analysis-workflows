#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "pindel parallel workflow"
requirements:
    - class: ScatterFeatureRequirement
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
    pindel_cat:
        scatter: interval_list
        run: pindel_cat.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: split_interval_list/split_interval_lists
            insert_size: insert_size
        out:
            [per_interval_pindel_out]
    cat_all:
        run: cat_all.cwl
        in:
            interval_pindel_outs: [pindel_cat/per_interval_pindel_out]
        out:
            [all_interval_pindel_out]
    grep:
        run: grep.cwl
        in: 
           pindel_output: cat_all/all_interval_pindel_out
        out:
           [pindel_head] 
    somaticfilter:
        run: somaticfilter.cwl
        in:
            reference: reference
            pindel_output_summary: grep/pindel_head
        out: 
            [vcf]
    bgzip:
        run: ../detect_variants/bgzip.cwl
        in: 
            file: somaticfilter/vcf
        out:
            [bgzipped_file]
    index:
        run: ../detect_variants/index.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
