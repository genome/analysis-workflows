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
    collect_hs_metrics_per_target_coverage:
        type: boolean?
    collect_hs_metrics_per_base_coverage:
        type: boolean?
    picard_metric_accumulation_level:
        type: string?
        default: ALL_READS
outputs:
    insert_size_metrics:
        type: File
        outputSource: collect_insert_size_metrics/insert_size_metrics
    alignment_summary_metrics:
        type: File
        outputSource: collect_alignment_summary_metrics/alignment_summary_metrics
    hs_metrics:
        type: File
        outputSource: collect_hs_metrics/hs_metrics
    per_target_coverage_metrics:
        type: File?
        outputSource: collect_hs_metrics/per_target_coverage_metrics
    per_base_coverage_metrics:
        type: File?
        outputSource: collect_hs_metrics/per_base_coverage_metrics
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
            [insert_size_metrics]
    collect_alignment_summary_metrics:
        run: collect_alignment_summary_metrics.cwl
        in:
            cram: cram
            reference: reference
            metric_accumulation_level: picard_metric_accumulation_level
        out:
            [alignment_summary_metrics]
    collect_hs_metrics:
        run: collect_hs_metrics.cwl
        in:
            cram: cram
            reference: reference
            metric_accumulation_level: picard_metric_accumulation_level
            bait_intervals: bait_intervals
            target_intervals: target_intervals
            per_target_coverage: collect_hs_metrics_per_target_coverage
            per_base_coverage: collect_hs_metrics_per_base_coverage
        out:
            [hs_metrics, per_target_coverage_metrics, per_base_coverage_metrics]
    samtools_flagstat:
        run: samtools_flagstat.cwl
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
