#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "panel of normals workflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type: File
        secondaryFiles: [".fai"]
    tumor_bams:
        type: File[]
    scatter_count:
        type: int
        default: 50
    dbsnp_vcf:
        type: File
        secondaryFiles: [.tbi]
    cosmic_vcf:
        type: File
        secondaryFiles: [.tbi]
    artifact_detection_mode:
        type: boolean
        default: true
    interval_list:
        type: File
    minimumN:
        type: int
        default: 2
outputs:
    pon_vcf:
        type: File
        outputSource: combine_variants/combined_vcf
        secondaryFiles: [.tbi]
steps:
    mutect:
        scatter: tumor_bam
        run: ../mutect/workflow.cwl
        in:
            reference: reference
            tumor_bam: tumor_bams
            dbsnp_vcf: dbsnp_vcf
            cosmic_vcf: cosmic_vcf
            artifact_detection_mode: artifact_detection_mode
            interval_list: interval_list
            scatter_count: scatter_count
        out:
            [merged_vcf]
    combine_variants:
        run: combine_variants.cwl
        in:
            reference: reference
            interval_list: interval_list
            minimumN: minimumN
            vcfs: [mutect/merged_vcf]
        out:
            [combined_vcf]
