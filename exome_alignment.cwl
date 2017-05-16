#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "exome alignment with qc"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type: File
        secondaryFiles: [.fai, .bwt, .sa, .ann, .amb, .pac, ^.dict, .alt]
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
    dbsnp:
        type: File
        secondaryFiles: [.tbi]
    bait_intervals:
        type: File
    target_intervals:
        type: File
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
outputs:
    cram:
        type: File
        outputSource: alignment/final_cram
    insert_size_metrics:
        type: File
        outputSource: qc/insert_size_metrics
    alignment_summary_metrics:
        type: File
        outputSource: qc/alignment_summary_metrics
    hs_metrics:
        type: File
        outputSource: qc/hs_metrics
    per_target_coverage_metrics:
        type: File
        outputSource: qc/per_target_coverage_metrics
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
        run: unaligned_bam_to_bqsr/workflow.cwl
        in:
            reference: reference
            bams: bams
            readgroups: readgroups
            mills: mills
            known_indels: known_indels
            dbsnp: dbsnp
        out: [final_cram]
    qc:
        run: qc/workflow_exome.cwl
        in:
            cram: alignment/final_cram
            reference: reference
            bait_intervals: bait_intervals
            target_intervals: target_intervals
            omni_vcf: omni_vcf
        out: [insert_size_metrics, alignment_summary_metrics, hs_metrics, per_target_coverage_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
