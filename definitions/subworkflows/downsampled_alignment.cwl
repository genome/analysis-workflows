#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "downsample unaligned BAM and align"
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement

inputs:
    bams:
        type: File[]
    readgroups:
        type: string[]
    reference:
        type: string
    dbsnp_vcf:
        type: File
        secondaryFiles: [.tbi]
    mills:
        type: File
        secondaryFiles: [.tbi]
    known_indels:
        type: File
        secondaryFiles: [.tbi]
    downsample_probability:
        type: float
outputs:
    cram:
        type: File
        outputSource: align_workflow/final_cram
        secondaryFiles: [.crai, ^.crai]
    mark_duplicates_metrics_file:
        type: File
        outputSource: align_workflow/mark_duplicates_metrics_file
steps:
    downsample:
        scatter: [bam]
        scatterMethod: dotproduct
        run: ../tools/downsample.cwl
        in:
            bam: bams
            probability: downsample_probability
            reference: reference
        out:
            [downsampled_bam]
    align_workflow:
        run: bam_to_bqsr.cwl
        in:
            bams: downsample/downsampled_bam
            readgroups: readgroups
            reference: reference
            dbsnp_vcf: dbsnp_vcf
            mills: mills
            known_indels: known_indels
        out:
            [final_cram,mark_duplicates_metrics_file]
