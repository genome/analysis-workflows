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
    bqsr_known_sites:
        type: File[]
        secondaryFiles: [.tbi]
        doc: "One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis."
    downsample_probability:
        type: float
outputs:
    bam:
        type: File
        outputSource: align_workflow/final_bam
        secondaryFiles: [.bai, ^.bai]
    mark_duplicates_metrics_file:
        type: File
        outputSource: align_workflow/mark_duplicates_metrics_file
steps:
    downsample:
        scatter: [sam]
        scatterMethod: dotproduct
        run: ../tools/downsample.cwl
        in:
            sam: bams
            probability: downsample_probability
            reference: reference
        out:
            [downsampled_sam]
    align_workflow:
        run: bam_to_bqsr.cwl
        in:
            bams: downsample/downsampled_sam
            readgroups: readgroups
            reference: reference
            bqsr_known_sites: bqsr_known_sites
        out:
            [final_bam,mark_duplicates_metrics_file]
