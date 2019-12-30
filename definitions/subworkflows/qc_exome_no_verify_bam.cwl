#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Exome QC workflow"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
    - class: SubworkflowFeatureRequirement
inputs:
    bam:
        type: File
        secondaryFiles: [^.bai]
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    bait_intervals:
        type: File
    target_intervals:
        type: File
    picard_metric_accumulation_level:
        type: string
        default: ALL_READS
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
    insert_size_metrics:
        type: File
        outputSource: collect_insert_size_metrics/insert_size_metrics
    insert_size_histogram:
        type: File
        outputSource: collect_insert_size_metrics/insert_size_histogram
    alignment_summary_metrics:
        type: File
        outputSource: collect_alignment_summary_metrics/alignment_summary_metrics
    hs_metrics:
        type: File
        outputSource: collect_roi_hs_metrics/hs_metrics
    per_target_coverage_metrics:
        type: File[]
        outputSource: collect_detailed_hs_metrics/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File[]
        outputSource: collect_detailed_hs_metrics/per_target_hs_metrics
    per_base_coverage_metrics:
        type: File[]
        outputSource: collect_detailed_hs_metrics/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File[]
        outputSource: collect_detailed_hs_metrics/per_base_hs_metrics
    summary_hs_metrics:
        type: File[]
        outputSource: collect_detailed_hs_metrics/summary_hs_metrics
    flagstats:
        type: File
        outputSource: samtools_flagstat/flagstats
steps:
    collect_insert_size_metrics:
        run: ../tools/collect_insert_size_metrics.cwl
        in:
            bam: bam
            reference: reference
            metric_accumulation_level: picard_metric_accumulation_level
        out:
            [insert_size_metrics, insert_size_histogram]
    collect_alignment_summary_metrics:
        run: ../tools/collect_alignment_summary_metrics.cwl
        in:
            bam: bam
            reference: reference
            metric_accumulation_level: picard_metric_accumulation_level
        out:
            [alignment_summary_metrics]
    collect_roi_hs_metrics:
        run: ../tools/collect_hs_metrics.cwl
        in:
            bam: bam
            reference: reference
            metric_accumulation_level:
                valueFrom: "ALL_READS"
            bait_intervals: bait_intervals
            target_intervals: target_intervals
            per_target_coverage:
                default: false
            per_base_coverage:
                default: false
            output_prefix:
                valueFrom: "roi"
            minimum_mapping_quality: minimum_mapping_quality
            minimum_base_quality: minimum_base_quality
        out:
            [hs_metrics]
    collect_detailed_hs_metrics:
        run: hs_metrics.cwl
        in:
            bam: bam
            minimum_mapping_quality: minimum_mapping_quality
            minimum_base_quality: minimum_base_quality
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            reference: reference
            summary_intervals: summary_intervals
        out:
            [per_base_coverage_metrics, per_base_hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, summary_hs_metrics]
    samtools_flagstat:
        run: ../tools/samtools_flagstat.cwl
        in:
            bam: bam
        out: [flagstats]
