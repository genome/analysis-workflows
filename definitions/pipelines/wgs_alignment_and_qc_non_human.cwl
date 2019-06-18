#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "wgs alignment with qc"
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
    intervals:
        type: File
    picard_metric_accumulation_level:
        type: string
    minimum_mapping_quality:
        type: int?
    minimum_base_quality:
        type: int?
    per_base_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    per_target_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    summary_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
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
    gc_bias_metrics:
        type: File
        outputSource: qc/gc_bias_metrics
    gc_bias_metrics_chart:
        type: File
        outputSource: qc/gc_bias_metrics_chart
    gc_bias_metrics_summary:
        type: File
        outputSource: qc/gc_bias_metrics_summary
    wgs_metrics:
        type: File
        outputSource: qc/wgs_metrics
    flagstats:
        type: File
        outputSource: qc/flagstats
    per_base_coverage_metrics:
        type: File[]
        outputSource: qc/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File[]
        outputSource: qc/per_base_hs_metrics
    per_target_coverage_metrics:
        type: File[]
        outputSource: qc/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File[]
        outputSource: qc/per_target_hs_metrics
    summary_hs_metrics:
        type: File[]
        outputSource: qc/summary_hs_metrics
    bamcoverage_bigwig:
        type: File
        outputSource: qc/bamcoverage_bigwig
steps:
    alignment:
        run: ../subworkflows/align_sort_markdup.cwl
        in:
            reference: reference
            bams: bams
            readgroups: readgroups
        out: [final_bam,mark_duplicates_metrics_file]
    qc:
        run: ../subworkflows/qc_wgs_no_verify_bam.cwl
        in:
            bam: alignment/final_bam
            reference: reference
            intervals: intervals
            picard_metric_accumulation_level: picard_metric_accumulation_level
            minimum_mapping_quality: minimum_mapping_quality
            minimum_base_quality: minimum_base_quality
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
        out: [insert_size_metrics, insert_size_histogram, alignment_summary_metrics, gc_bias_metrics, gc_bias_metrics_chart, gc_bias_metrics_summary, wgs_metrics, flagstats, per_base_coverage_metrics, per_base_hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, summary_hs_metrics, bamcoverage_bigwig]
