#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "umi molecular alignment workflow"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
inputs:
    bam:
        type: File[]
    sample_name:
        type: string
    read_structure:
        type: string[]
    reference:
        type: string
    target_intervals:
        type: File
    bait_intervals:
        type: File
    per_target_intervals:
        type: File
    per_target_bait_intervals:
        type: File
    per_base_intervals:
        type: File
    per_base_bait_intervals:
        type: File
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
    picard_metric_accumulation_level:
        type: string
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
    aligned_cram:
        type: File
        secondaryFiles: [.crai, ^.crai]
        outputSource: alignment/aligned_cram
    adapter_histogram:
        type: File[]
        outputSource: alignment/adapter_histogram
    duplex_seq_metrics:
        type: File[]
        outputSource: alignment/duplex_seq_metrics
    insert_size_histogram:
        type: File
        outputSource: qc/insert_size_histogram
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
    verify_bam_id_metrics:
        type: File
        outputSource: qc/verify_bam_id_metrics
    verify_bam_id_depth:
        type: File
        outputSource: qc/verify_bam_id_depth
steps:
    alignment:
        run: molecular_alignment.cwl
        in:
            bam: bam
            sample_name: sample_name
            read_structure: read_structure
            reference: reference
            target_intervals: target_intervals
        out:
            [aligned_cram, adapter_histogram, duplex_seq_metrics]
    qc:
        run: qc_exome.cwl
        in:
            cram: alignment/aligned_cram
            reference: reference
            bait_intervals: bait_intervals
            target_intervals: target_intervals
            omni_vcf: omni_vcf
            picard_metric_accumulation_level: picard_metric_accumulation_level
            minimum_mapping_quality: minimum_mapping_quality
            minimum_base_quality: minimum_base_quality
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
        out:
            [insert_size_metrics, insert_size_histogram, alignment_summary_metrics, hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, per_base_coverage_metrics, per_base_hs_metrics, summary_hs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
