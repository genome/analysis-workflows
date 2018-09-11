#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Exome QC workflow"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    cram:
        type: File
        secondaryFiles: [^.crai]
    reference:
        type: string
    bait_intervals:
        type: File
    target_intervals:
        type: File
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
    per_target_intervals:
        type: File
    per_target_bait_intervals:
        type: File
    per_base_intervals:
        type: File
    per_base_bait_intervals:
        type: File
    picard_metric_accumulation_level:
        type: string?
        default: ALL_READS
    minimum_mapping_quality:
        type: int?
    minimum_base_quality:
        type: int?
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
        type: File?
        outputSource: collect_per_target_hs_metrics/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File?
        outputSource: collect_per_target_hs_metrics/hs_metrics
    per_base_coverage_metrics:
        type: File?
        outputSource: collect_per_base_hs_metrics/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File?
        outputSource: collect_per_base_hs_metrics/hs_metrics
    flagstats:
        type: File
        outputSource: samtools_flagstat/flagstats
    verify_bam_id_metrics:
        type: File
        outputSource: verify_bam_id/verify_bam_id_metrics
    verify_bam_id_depth:
        type: File
        outputSource: verify_bam_id/verify_bam_id_depth
steps:
    collect_insert_size_metrics:
        run: collect_insert_size_metrics.cwl
        in:
            cram: cram
            reference: reference
            metric_accumulation_level: picard_metric_accumulation_level
        out:
            [insert_size_metrics, insert_size_histogram]
    collect_alignment_summary_metrics:
        run: collect_alignment_summary_metrics.cwl
        in:
            cram: cram
            reference: reference
            metric_accumulation_level: picard_metric_accumulation_level
        out:
            [alignment_summary_metrics]
    collect_roi_hs_metrics:
        run: collect_hs_metrics.cwl
        in:
            cram: cram
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
    collect_per_base_hs_metrics:
        run: collect_hs_metrics.cwl
        in:
            cram: cram
            reference: reference
            metric_accumulation_level:
                valueFrom: "ALL_READS"
            bait_intervals: per_base_bait_intervals
            target_intervals: per_base_intervals
            per_target_coverage:
                default: false
            per_base_coverage:
                default: true
            output_prefix:
                valueFrom: "base"
            minimum_mapping_quality: minimum_mapping_quality
            minimum_base_quality: minimum_base_quality
        out:
            [hs_metrics, per_base_coverage_metrics]
    collect_per_target_hs_metrics:
        run: collect_hs_metrics.cwl
        in:
            cram: cram
            reference: reference
            metric_accumulation_level:
                valueFrom: "ALL_READS"
            bait_intervals: per_target_bait_intervals
            target_intervals: per_target_intervals
            per_target_coverage:
                default: true
            per_base_coverage:
                default: false
            output_prefix:
                valueFrom: "target"
            minimum_mapping_quality: minimum_mapping_quality
            minimum_base_quality: minimum_base_quality
        out:
            [hs_metrics, per_target_coverage_metrics]
    samtools_flagstat:
        run: ../definitions/tools/samtools_flagstat.cwl
        in:
            cram: cram
        out: [flagstats]
    select_variants:
        run: ../detect_variants/select_variants.cwl
        in:
            reference: reference
            vcf: omni_vcf
            interval_list: target_intervals
        out:
            [filtered_vcf]
    cram_to_bam:
        run: ../cram_to_bam/workflow.cwl
        in:
          cram: cram
          reference: reference
        out:
          [bam]
    verify_bam_id:
        run: verify_bam_id.cwl
        in:
            bam: cram_to_bam/bam
            vcf: select_variants/filtered_vcf
        out:
            [verify_bam_id_metrics, verify_bam_id_depth]
