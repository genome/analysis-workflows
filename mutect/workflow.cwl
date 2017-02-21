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
        secondaryFiles: [^.bai]
    normal_bam:
        type: File?
        secondaryFiles: [^.bai]
    interval_list:
        type: File
    scatter_count:
        type: int
        default: 50
    dbsnp_vcf:
        type: File?
        secondaryFiles: [.tbi]
    cosmic_vcf:
        type: File?
        secondaryFiles: [.tbi]
    artifact_detection_mode:
        type: boolean?
outputs:
    merged_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [.tbi]
steps:
    split_interval_list:
        run: ../detect_variants/split_interval_list.cwl
        in:
            interval_list: interval_list
            scatter_count: scatter_count
        out: [split_interval_lists]
    mutect:
        scatter: interval_list
        run: mutect.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: split_interval_list/split_interval_lists
            dbsnp_vcf: dbsnp_vcf
            cosmic_vcf: cosmic_vcf
            artifact_detection_mode: artifact_detection_mode
        out:
            [vcf]
    merge:
        run: ../detect_variants/merge.cwl
        in:
            vcfs: [mutect/vcf]
        out:
            [merged_vcf]
    index:
        run: ../detect_variants/index.cwl
        in:
            vcf: [merge/merged_vcf]
        out:
            [indexed_vcf]

