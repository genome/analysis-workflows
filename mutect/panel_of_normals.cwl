#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "mutect panel-of-normals workflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type: File
        secondaryFiles: [".fai", "^.dict"]
    normal_crams:
        type: File[]
        secondaryFiles: [^.crai]
    interval_list:
        type: File
    scatter_count:
        type: int
    dbsnp_vcf:
        type: File?
        secondaryFiles: [.tbi]
    cosmic_vcf:
        type: File?
        secondaryFiles: [.tbi]
    artifact_detection_mode:
        type: boolean
        default: true
outputs:
    pon_vcf:
        type: File
        outputSource: combine/combined_vcf
        secondaryFiles: [.tbi]
steps:
    mutect:
        scatter: tumor_cram
        run: workflow.cwl
        in:
            reference: reference
            tumor_cram: normal_crams
            interval_list: interval_list
            dbsnp_vcf: dbsnp_vcf
            cosmic_vcf: cosmic_vcf
            artifact_detection_mode: artifact_detection_mode
            scatter_count: scatter_count
        out:
            [merged_vcf]
    combine:
        run: combine_pon.cwl
        in:
            reference: reference
            interval_list: interval_list
            normal_vcfs: mutect/merged_vcf
        out:
            [combined_vcf]

