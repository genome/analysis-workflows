#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "mutect parallel workflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type: File
        secondaryFiles: [".fai", "^.dict"]
    tumor_bam:
        type: File
        secondaryFiles: .bai
    normal_bam:
        type: File
        secondaryFiles: .bai
    interval_list:
        type: File[]
outputs:
    merged_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: .tbi
steps:
    mutect_and_index:
        scatter: interval_list
        run: mutect_and_index.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: interval_list
        out:
            [vcf]
    merge:
        run: ../detect_variants/merge.cwl
        in:
            vcfs: [mutect_and_index/vcf]
        out:
            [merged_vcf]
    index:
        run: ../detect_variants/index.cwl
        in:
            vcf: [merge/merged_vcf]
        out:
            [indexed_vcf]

