#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Chipseq alignment with qc and creating homer tag directory"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
          - $import: ../types/sequence_data.yml
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
    final_name:
        type: string?
    chipseq_sequence:
        type: ../types/sequence_data.yml#sequence_data[]
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
    intervals:
        type: File
    per_base_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
        default: []
    per_target_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
        default: []
    summary_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
        default: []
    picard_metric_accumulation_level:
        type: string
    bqsr_known_sites:
        type: File[]
        secondaryFiles: [.tbi]
        doc: "One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis."
    bqsr_intervals:
        type: string[]
    minimum_mapping_quality:
        type: int?
    minimum_base_quality:
        type: int?
outputs:
    bam:
        type: File
        outputSource: alignment/final_bam
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
    per_base_coverage_metrics:
        type: File[]
        outputSource: qc/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File[]
        outputSource: qc/per_base_hs_metrics
    per_target_coverage_metrics:
        type: File[]
        outputSource: qc/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File[]
        outputSource: qc/per_target_hs_metrics
    summary_hs_metrics:
        type: File[]
        outputSource: qc/summary_hs_metrics
    tag_directory:
        type: Directory
        outputSource: homer_tag_directory/tag_directory

steps:
    alignment:
        run: ../subworkflows/sequence_to_bqsr.cwl
        in:
            reference: reference
            unaligned: chipseq_sequence
            bqsr_known_sites: bqsr_known_sites
            bqsr_intervals: bqsr_intervals
            final_name: final_name
        out: [final_bam,mark_duplicates_metrics_file]

    bam_to_sam:
        run: ../tools/bam_to_sam.cwl
        in:
            bam: alignment/final_bam
        out: [final_sam]

    homer_tag_directory:
        run: ../tools/homer_tag_directory.cwl
        in:
            sam: bam_to_sam/final_sam
        out: [tag_directory]

    qc:
        run: ../subworkflows/qc_wgs.cwl
        in:
            bam: alignment/final_bam
            reference: reference
            omni_vcf: omni_vcf
            intervals: intervals
            picard_metric_accumulation_level: picard_metric_accumulation_level
            minimum_mapping_quality: minimum_mapping_quality
            minimum_base_quality: minimum_base_quality
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
        out: [insert_size_metrics, insert_size_histogram, alignment_summary_metrics, gc_bias_metrics, gc_bias_metrics_chart, gc_bias_metrics_summary, wgs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth, per_base_coverage_metrics, per_base_hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, summary_hs_metrics]
