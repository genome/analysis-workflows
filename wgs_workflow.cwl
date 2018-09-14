#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "wgs alignment and germline variant detection"
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
    qc_intervals:
        type: File
    picard_metric_accumulation_level:
        type: string
    variant_detection_intervals:
        type: File
    vep_cache_dir:
        type: string?
    synonyms_file:
        type: File?
    sample_name:
        type: string
    docm_vcf:
        type: File
    custom_gnomad_vcf:
        type: File?
        secondaryFiles: [.tbi]
    readcount_minimum_mapping_quality:
        type: int?
    readcount_minimum_base_quality:
        type: int?
outputs:
    cram:
        type: File
        outputSource: alignment_and_qc/cram
    mark_duplicates_metrics:
        type: File
        outputSource: alignment_and_qc/mark_duplicates_metrics
    insert_size_metrics:
        type: File
        outputSource: alignment_and_qc/insert_size_metrics
    insert_size_histogram:
        type: File
        outputSource: alignment_and_qc/insert_size_histogram
    alignment_summary_metrics:
        type: File
        outputSource: alignment_and_qc/alignment_summary_metrics
    gc_bias_metrics:
        type: File
        outputSource: alignment_and_qc/gc_bias_metrics
    gc_bias_metrics_chart:
        type: File
        outputSource: alignment_and_qc/gc_bias_metrics_chart
    gc_bias_metrics_summary:
        type: File
        outputSource: alignment_and_qc/gc_bias_metrics_summary
    wgs_metrics:
        type: File
        outputSource: alignment_and_qc/wgs_metrics
    flagstats:
        type: File
        outputSource: alignment_and_qc/flagstats
    verify_bam_id_metrics:
        type: File
        outputSource: alignment_and_qc/verify_bam_id_metrics
    verify_bam_id_depth:
        type: File
        outputSource: alignment_and_qc/verify_bam_id_depth
    varscan_vcf:
        type: File
        outputSource: detect_variants/varscan_vcf
        secondaryFiles: [.tbi]
    docm_gatk_vcf:
        type: File
        outputSource: detect_variants/docm_gatk_vcf
    annotated_vcf:
        type: File
        outputSource: detect_variants/annotated_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: detect_variants/final_vcf
        secondaryFiles: [.tbi]
    final_tsv:
        type: File
        outputSource: detect_variants/final_tsv
    vep_summary:
        type: File
        outputSource: detect_variants/vep_summary
    tumor_bam_readcount_tsv:
        type: File
        outputSource: detect_variants/tumor_bam_readcount_tsv
steps:
    alignment_and_qc:
        run: definitions/pipelines/wgs_alignment.cwl
        in:
            reference: reference
            bams: bams
            readgroups: readgroups
            mills: mills
            known_indels: known_indels
            dbsnp_vcf: dbsnp_vcf
            omni_vcf: omni_vcf
            intervals: qc_intervals
            picard_metric_accumulation_level: picard_metric_accumulation_level
        out:
            [cram, mark_duplicates_metrics, insert_size_metrics, insert_size_histogram, alignment_summary_metrics, gc_bias_metrics, gc_bias_metrics_chart, gc_bias_metrics_summary, wgs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
    detect_variants:
        run: detect_variants/tumor_only_detect_variants.cwl
        in:
            reference: reference
            cram: alignment_and_qc/cram
            interval_list: variant_detection_intervals
            #varscan_strand_filter:
            #varscan_min_coverage:
            #varscan_min_var_freq:
            #varscan_p_value:
            #varscan_min_reads:
            #maximum_population_allele_frequency:
            vep_cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            #variants_to_table_fields:
            #variants_to_table_genotype_fields:
            #vep_to_table_fields:
            sample_name: sample_name
            docm_vcf: docm_vcf
            custom_gnomad_vcf: custom_gnomad_vcf
            readcount_minimum_mapping_quality: readcount_minimum_mapping_quality
            readcount_minimum_base_quality: readcount_minimum_base_quality
        out:
            [varscan_vcf, docm_gatk_vcf, annotated_vcf, final_vcf, final_tsv, vep_summary, tumor_bam_readcount_tsv]
