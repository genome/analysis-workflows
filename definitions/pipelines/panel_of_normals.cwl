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
outputs:
    pon_vcf:
        type: File
        outputSource: combine/combined_vcf
        secondaryFiles: [.tbi]
steps:
    mutect:
        scatter: tumor_cram
        run: ../subworkflows/mutect.cwl
        in:
            reference: reference
            tumor_cram: normal_crams
            normal_cram:
                default: null
            interval_list: interval_list
            scatter_count: scatter_count
            dbsnp_vcf: dbsnp_vcf
            cosmic_vcf: cosmic_vcf
            artifact_detection_mode:
                default: true
            panel_of_normals_vcf:
                default: null
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

