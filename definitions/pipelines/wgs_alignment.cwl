#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "wgs alignment with qc"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference: string
    bams:
        type: File[]
    readgroups:
        type: string[]
    mills:
        type: File
        secondaryFiles: [.tbi]
    known_indels:
        type: File
        secondaryFiles: [.tbi]
    dbsnp_vcf:
        type: File
        secondaryFiles: [.tbi]
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
    intervals:
        type: File
    picard_metric_accumulation_level:
        type: string
outputs:
    cram:
        type: File
        outputSource: alignment/final_cram
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
    verify_bam_id_metrics:
        type: File
        outputSource: qc/verify_bam_id_metrics
    verify_bam_id_depth:
        type: File
        outputSource: qc/verify_bam_id_depth
steps:
    alignment:
        run: ../subworkflows/bam_to_bqsr_workflow.cwl
        in:
            reference: reference
            bams: bams
            readgroups: readgroups
            mills: mills
            known_indels: known_indels
            dbsnp_vcf: dbsnp_vcf
        out: [final_cram,mark_duplicates_metrics_file]
    qc:
        run: ../subworkflows/qc_wgs.cwl
        in:
            cram: alignment/final_cram
            reference: reference
            omni_vcf: omni_vcf
            intervals: intervals
            picard_metric_accumulation_level: picard_metric_accumulation_level
        out: [insert_size_metrics, insert_size_histogram, alignment_summary_metrics, gc_bias_metrics, gc_bias_metrics_chart, gc_bias_metrics_summary, wgs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
