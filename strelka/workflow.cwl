#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "strelka workflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement
inputs:
    tumor_bam:
        type: File
        secondaryFiles: [^.bai]
    normal_bam:
        type: File
        secondaryFiles: [^.bai]
    reference:
        type: File
        secondaryFiles: [.fai]
    interval_list:
        type: File
    exome_mode:
        type: boolean
outputs:
    merged_vcf:
        type: File
        outputSource: index_filtered/indexed_vcf
        secondaryFiles: [.tbi]
steps:
    strelka:
        run: strelka.cwl
        in:
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            reference: reference
            exome_mode: exome_mode
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
            vcfs: [process/processed_vcf]
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
    index_filtered:
        run: ../detect_variants/index.cwl
        in:
            vcf: region_filter/filtered_vcf
        out:
            [indexed_vcf]



