#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "exome alignment and germline variant detection, with optitype for HLA typing"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
          - $import: ../types/sequence_data.yml
          - $import: ../types/vep_custom_annotation.yml
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
    sequence:
        type: ../types/sequence_data.yml#sequence_data[]
    mills:
        type: File
        secondaryFiles: [.tbi]
    known_indels:
        type: File
        secondaryFiles: [.tbi]
    dbsnp_vcf:
        type: File
        secondaryFiles: [.tbi]
    bqsr_intervals:
        type: string[]?
    bait_intervals:
        type: File
    target_intervals:
        type: File
    per_base_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    per_target_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    summary_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
    picard_metric_accumulation_level:
        type: string
    emit_reference_confidence:
        type:
            type: enum
            symbols: ['NONE', 'BP_RESOLUTION', 'GVCF']
    gvcf_gq_bands:
        type: string[]
    intervals:
        type:
            type: array
            items:
                type: array
                items: string
    vep_cache_dir:
        type: string
    vep_ensembl_assembly:
        type: string
        doc: "genome assembly to use in vep. Examples: GRCh38 or GRCm38"
    vep_ensembl_version:
        type: string
        doc: "ensembl version - Must be present in the cache directory. Example: 95"
    vep_ensembl_species:
        type: string
        doc: "ensembl species - Must be present in the cache directory. Examples: homo_sapiens or mus_musculus"
    synonyms_file:
        type: File?
    annotate_coding_only:
        type: boolean?
    vep_custom_annotations:
        type: ../types/vep_custom_annotation.yml#vep_custom_annotation[]
        doc: "custom type, check types directory for input format"
    qc_minimum_mapping_quality:
        type: int?
    qc_minimum_base_quality:
        type: int?
    optitype_name:
        type: string?
outputs:
    cram:
        type: File
        outputSource: germline_exome/cram
    mark_duplicates_metrics:
        type: File
        outputSource: germline_exome/mark_duplicates_metrics
    insert_size_metrics:
        type: File
        outputSource: germline_exome/insert_size_metrics
    insert_size_histogram:
        type: File
        outputSource: germline_exome/insert_size_histogram
    alignment_summary_metrics:
        type: File
        outputSource: germline_exome/alignment_summary_metrics
    hs_metrics:
        type: File
        outputSource: germline_exome/hs_metrics
    per_target_coverage_metrics:
        type: File[]
        outputSource: germline_exome/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File[]
        outputSource: germline_exome/per_target_hs_metrics
    per_base_coverage_metrics:
        type: File[]
        outputSource: germline_exome/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File[]
        outputSource: germline_exome/per_base_hs_metrics
    summary_hs_metrics:
        type: File[]
        outputSource: germline_exome/summary_hs_metrics
    flagstats:
        type: File
        outputSource: germline_exome/flagstats
    verify_bam_id_metrics:
        type: File
        outputSource: germline_exome/verify_bam_id_metrics
    verify_bam_id_depth:
        type: File
        outputSource: germline_exome/verify_bam_id_depth
    gvcf:
        type: File[]
        outputSource: germline_exome/gvcf
    final_vcf:
        type: File
        outputSource: germline_exome/final_vcf
        secondaryFiles: [.tbi]
    filtered_vcf:
        type: File
        outputSource: germline_exome/filtered_vcf
        secondaryFiles: [.tbi]
    vep_summary:
        type: File
        outputSource: germline_exome/vep_summary
    optitype_tsv:
        type: File
        outputSource: optitype/optitype_tsv
    optitype_plot:
        type: File
        outputSource: optitype/optitype_plot
steps:
    germline_exome:
        run: germline_exome.cwl
        in:
            reference: reference
            sequence: sequence
            mills: mills
            known_indels: known_indels
            dbsnp_vcf: dbsnp_vcf
            bqsr_intervals: bqsr_intervals
            bait_intervals: bait_intervals
            target_intervals: target_intervals
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
            omni_vcf: omni_vcf
            picard_metric_accumulation_level: picard_metric_accumulation_level
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            vep_cache_dir: vep_cache_dir
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            vep_custom_annotations: vep_custom_annotations
            qc_minimum_mapping_quality: qc_minimum_mapping_quality
            qc_minimum_base_quality: qc_minimum_base_quality
        out:
            [cram, mark_duplicates_metrics, insert_size_metrics, insert_size_histogram, alignment_summary_metrics, hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, per_base_coverage_metrics, per_base_hs_metrics, summary_hs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth, gvcf, final_vcf, filtered_vcf, vep_summary]
    optitype:
        run: ../tools/optitype_dna.cwl
        in:
            optitype_name: optitype_name
            cram: germline_exome/cram
        out:
            [optitype_tsv, optitype_plot]
