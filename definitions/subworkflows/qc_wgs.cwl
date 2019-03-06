#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "WGS QC workflow"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
    - class: SubworkflowFeatureRequirement
inputs:
    cram:
        type: File
        secondaryFiles: [^.crai]
    reference:
        type: string
    intervals:
        type: File
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
    picard_metric_accumulation_level:
        type: string?
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
    gc_bias_metrics:
        type: File
        outputSource: collect_gc_bias_metrics/gc_bias_metrics
    gc_bias_metrics_chart:
        type: File
        outputSource: collect_gc_bias_metrics/gc_bias_metrics_chart
    gc_bias_metrics_summary:
        type: File
        outputSource: collect_gc_bias_metrics/gc_bias_metrics_summary
    wgs_metrics:
        type: File
        outputSource: collect_wgs_metrics/wgs_metrics
    flagstats:
        type: File
        outputSource: samtools_flagstat/flagstats
    verify_bam_id_metrics:
        type: File
        outputSource: verify_bam_id/verify_bam_id_metrics
    verify_bam_id_depth:
        type: File
        outputSource: verify_bam_id/verify_bam_id_depth
    per_base_coverage_metrics:
        type: File[]
        outputSource: collect_hs_metrics/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File[]
        outputSource: collect_hs_metrics/per_base_hs_metrics
    per_target_coverage_metrics:
        type: File[]
        outputSource: collect_hs_metrics/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File[]
        outputSource: collect_hs_metrics/per_target_hs_metrics
    summary_hs_metrics:
        type: File[]
        outputSource: collect_hs_metrics/summary_hs_metrics
    bamcoverage_bigwig:
        type: File
        outputSource: deeptools_bamcoverage/outfile
steps:
    collect_insert_size_metrics:
        run: ../tools/collect_insert_size_metrics.cwl
        in:
            cram: cram
            reference: reference
            metric_accumulation_level: picard_metric_accumulation_level
        out:
            [insert_size_metrics, insert_size_histogram]
    collect_alignment_summary_metrics:
        run: ../tools/collect_alignment_summary_metrics.cwl
        in:
            cram: cram
            reference: reference
            metric_accumulation_level: picard_metric_accumulation_level
        out:
            [alignment_summary_metrics]
    collect_gc_bias_metrics:
        run: ../tools/collect_gc_bias_metrics.cwl
        in:
            cram: cram
            reference: reference
            metric_accumulation_level: picard_metric_accumulation_level
        out:
            [gc_bias_metrics, gc_bias_metrics_chart, gc_bias_metrics_summary]
    collect_wgs_metrics:
        run: ../tools/collect_wgs_metrics.cwl
        in:
            cram: cram
            reference: reference
            intervals: intervals
        out:
            [wgs_metrics]
    samtools_flagstat:
        run: ../tools/samtools_flagstat.cwl
        in:
            cram: cram
        out: [flagstats]
    cram_to_bam:
        run: cram_to_bam_and_index.cwl
        in:
          cram: cram
          reference: reference
        out:
          [bam]
    verify_bam_id:
        run: ../tools/verify_bam_id.cwl
        in:
            bam: cram_to_bam/bam
            vcf: omni_vcf
        out:
            [verify_bam_id_metrics, verify_bam_id_depth]
    collect_hs_metrics:
        run: hs_metrics.cwl
        in:
            cram: cram
            minimum_mapping_quality: minimum_mapping_quality
            minimum_base_quality: minimum_base_quality
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            reference: reference
            summary_intervals: summary_intervals
        out:
            [per_base_coverage_metrics, per_base_hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, summary_hs_metrics]
    deeptools_bamcoverage:
        run: ../tools/deeptools_bamcoverage.cwl
        in:
            bam: cram
            min_mapping_quality: minimum_mapping_quality
        out:
            [outfile]
