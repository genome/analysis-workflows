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
        type: string
    normal_bams:
        type: File[]
        secondaryFiles: [^.bai]
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
    normal_sample_name:
        type: string
    tumor_sample_name:
        type: string
outputs:
    pon_vcf:
        type: File
        outputSource: combine/combined_vcf
        secondaryFiles: [.tbi]
steps:
    mutect:
        scatter: tumor_bam
        run: ../subworkflows/mutect.cwl
        in:
            reference: reference
            tumor_bam: normal_bams
            normal_bam:
                default: null
            interval_list: interval_list
            scatter_count: scatter_count
            dbsnp_vcf: dbsnp_vcf
            cosmic_vcf: cosmic_vcf
            artifact_detection_mode:
                default: true
            panel_of_normals_vcf:
                default: null
            tumor_sample_name: tumor_sample_name
        out:
            [unfiltered_vcf]
    combine:
        run: ../tools/pon_combine_variants.cwl
        in:
            reference: reference
            interval_list: interval_list
            normal_vcfs: mutect/unfiltered_vcf
        out:
            [combined_vcf]

