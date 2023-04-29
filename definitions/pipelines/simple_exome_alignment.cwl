#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "simple xome alignment with qc"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
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
    cram:
        type: File
        outputSource: index_cram/indexed_cram
    mark_duplicates_metrics:
        type: File
        outputSource: mark_duplicates_and_sort/metrics_file
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
    align:
        scatter: [bam, readgroup]
        scatterMethod: dotproduct
        run: ../subworkflows/align.cwl
        in:
            bam: bams
            readgroup: readgroups
            reference: reference
        out:
            [tagged_bam]
    merge:
        run: ../tools/merge_bams_samtools.cwl
        in:
            bams: align/tagged_bam
            name: final_name
        out:
            [merged_bam]
    name_sort:
        run: ../tools/name_sort.cwl
        in:
            bam: merge/merged_bam
        out:
            [name_sorted_bam]
    mark_duplicates_and_sort:
        run: ../tools/mark_duplicates_and_sort.cwl
        in:
            bam: name_sort/name_sorted_bam
        out:
            [sorted_bam, metrics_file]
    index_bam:
        run: ../tools/index_bam.cwl
        in:
            bam: mark_duplicates_and_sort/sorted_bam
        out:
            [indexed_bam]
    qc:
        run: ../subworkflows/qc_exome_no_vbi.cwl
        in:
            bam: index_bam/indexed_bam
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
    bam_to_cram:
        run: ../tools/bam_to_cram.cwl
        in:
            bam: index_bam/indexed_bam
            reference: reference
        out:
            [cram]
    index_cram:
         run: ../tools/index_cram.cwl
         in:
            cram: bam_to_cram/cram
         out:
            [indexed_cram]
