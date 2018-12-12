#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Subworkflow to allow calling cnvkit with cram instead of bam files"
inputs:
    normal_cram:
        type: File
    tumor_cram:
        type: File
    reference:
        type: string
    bait_intervals:
        type: File
    access:
        type: File?
    method:
        type: string?
    diagram:
        type: boolean?
    scatter_plot:
        type: boolean?
    drop_low_coverage:
        type: boolean?
outputs:
    intervals_antitarget:
        type: File?
        outputSource: run_cnvkit/intervals_antitarget
    intervals_target:
        type: File?
        outputSource: run_cnvkit/intervals_target
    normal_antitarget_coverage:
        type: File
        outputSource: run_cnvkit/normal_antitarget_coverage
    normal_target_coverage:
        type: File
        outputSource: run_cnvkit/normal_target_coverage
    reference_coverage:
        type: File?
        outputSource: run_cnvkit/reference_coverage
    cn_diagram:
        type: File?
        outputSource: run_cnvkit/cn_diagram
    cn_scatter_plot:
        type: File?
        outputSource: run_cnvkit/cn_scatter_plot
    tumor_antitarget_coverage:
        type: File
        outputSource: run_cnvkit/tumor_antitarget_coverage
    tumor_target_coverage:
        type: File
        outputSource: run_cnvkit/tumor_target_coverage
    tumor_bin_level_ratios:
        type: File
        outputSource: run_cnvkit/tumor_bin_level_ratios
    tumor_segmented_ratios:
        type: File
        outputSource: run_cnvkit/tumor_segmented_ratios
steps:
    normal_cram_to_bam:
        run: ../tools/cram_to_bam.cwl
        in:
            cram: normal_cram
            reference: reference
        out:
            [bam]
    tumor_cram_to_bam:
        run: ../tools/cram_to_bam.cwl
        in:
            cram: tumor_cram
            reference: reference
        out:
            [bam]
    run_cnvkit:
        run: ../tools/cnvkit_batch.cwl
        in:
            tumor_bam: tumor_cram_to_bam/bam
            normal_bam: normal_cram_to_bam/bam
            bait_intervals: bait_intervals
            reference: reference
            access: access
            method: method
            diagram: diagram
            scatter_plot: scatter_plot
            drop_low_coverage: drop_low_coverage
        out:
            [intervals_antitarget, intervals_target, normal_antitarget_coverage, normal_target_coverage, reference_coverage, cn_diagram, cn_scatter_plot, tumor_antitarget_coverage, tumor_target_coverage, tumor_bin_level_ratios, tumor_segmented_ratios]
