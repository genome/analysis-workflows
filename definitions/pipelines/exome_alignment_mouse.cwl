#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "exome alignment with qc, no bqsr, no verify_bam_id"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
    - class: SubworkflowFeatureRequirement
inputs:
    reference: string
    bams:
        type: File[]
    readgroups:
        type: string[]
    bait_intervals:
        type: File
    final_name:
        type: string?
    target_intervals:
        type: File
    per_base_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    per_target_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    summary_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    picard_metric_accumulation_level:
        type: string
    qc_minimum_mapping_quality:
        type: int?
    qc_minimum_base_quality:
        type: int?
outputs:
    bam:
        type: File
        outputSource: alignment/final_bam
    mark_duplicates_metrics:
        type: File
        outputSource: alignment/mark_duplicates_metrics_file
    insert_size_metrics:
        type: File
        outputSource: qc/insert_size_metrics
    insert_size_histogram:
        type: File
        outputSource: qc/insert_size_histogram
    alignment_summary_metrics:
        type: File
        outputSource: qc/alignment_summary_metrics
    hs_metrics:
        type: File
        outputSource: qc/hs_metrics
    per_target_coverage_metrics:
        type: File[]
        outputSource: qc/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File[]
        outputSource: qc/per_target_hs_metrics
    per_base_coverage_metrics:
        type: File[]
        outputSource: qc/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File[]
        outputSource: qc/per_base_hs_metrics
    summary_hs_metrics:
        type: File[]
        outputSource: qc/summary_hs_metrics
    flagstats:
        type: File
        outputSource: qc/flagstats
steps:
    alignment:
        run: ../subworkflows/align_sort_markdup.cwl
        in:
            reference: reference
            bams: bams
            readgroups: readgroups
            final_name: final_name
        out: [final_bam,mark_duplicates_metrics_file]
    qc:
        run: ../subworkflows/qc_exome_no_verify_bam.cwl
        in:
            bam: alignment/final_bam
            reference: reference
            bait_intervals: bait_intervals
            target_intervals: target_intervals
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
            picard_metric_accumulation_level: picard_metric_accumulation_level
            minimum_mapping_quality: qc_minimum_mapping_quality
            minimum_base_quality: qc_minimum_base_quality
        out: [insert_size_metrics, insert_size_histogram, alignment_summary_metrics, hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, per_base_coverage_metrics, per_base_hs_metrics, summary_hs_metrics, flagstats]

