#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "WGS QC workflow"
requirements:
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
