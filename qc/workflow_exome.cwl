#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Exome QC workflow"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    bam:
        type: File
        secondaryFiles: [.bai]
    reference:
        type: File
        secondaryFiles: [.fai]
    bait_intervals:
        type: File
    target_intervals:
        type: File
    roi_intervals:
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
    hs_metrics:
        type: File
        outputSource: collect_hs_metrics/hs_metrics
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
    collect_hs_metrics:
        run: collect_hs_metrics.cwl
        in:
            bam: bam
            reference: reference
            bait_intervals: bait_intervals
            target_intervals: target_intervals
        out:
            [hs_metrics]
    samtools_flagstat:
        run: samtools_flagstat.cwl
        in:
            bam: bam
        out: [flagstats]
    select_variants:
        run: ../detect_variants/select_variants.cwl
        in:
            reference: reference
            vcf: omni_vcf
            interval_list: roi_intervals
        out:
            [filtered_vcf]
    verify_bam_id:
        run: verify_bam_id.cwl
        in:
            bam: bam
            vcf: select_variants/filtered_vcf
        out:
            [verify_bam_id_metrics, verify_bam_id_depth]
