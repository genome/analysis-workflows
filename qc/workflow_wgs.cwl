#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "WGS QC workflow"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    bam:
        type: File
        secondaryFiles: [^.bai]
    reference:
        type: File
        secondaryFiles: [.fai]
    intervals:
        type: File
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
outputs:
    insert_size_metrics:
        type: File
        outputSource: collect_insert_size_metrics/insert_size_metrics
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
        run: collect_insert_size_metrics.cwl
        in:
            bam: bam
        out:
            [insert_size_metrics]
    collect_alignment_summary_metrics:
        run: collect_alignment_summary_metrics.cwl
        in:
            bam: bam
            reference: reference
        out:
            [alignment_summary_metrics]
    collect_gc_bias_metrics:
        run: collect_gc_bias_metrics.cwl
        in:
            bam: bam
            reference: reference
        out:
            [gc_bias_metrics, gc_bias_metrics_chart, gc_bias_metrics_summary]
    collect_wgs_metrics:
        run: collect_wgs_metrics.cwl
        in:
            bam: bam
            reference: reference
            intervals: intervals
        out:
            [wgs_metrics]
    samtools_flagstat:
        run: samtools_flagstat.cwl
        in:
            bam: bam
        out: [flagstats]
    verify_bam_id:
        run: verify_bam_id.cwl
        in:
            bam: bam
            vcf: omni_vcf
        out:
            [verify_bam_id_metrics, verify_bam_id_depth]
