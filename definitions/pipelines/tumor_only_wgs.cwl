#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "wgs alignment and tumor-only variant detection"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
          - $import: ../types/sequence_data.yml
          - $import: ../types/trimming_options.yml
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
    trimming:
        type:
            - ../types/trimming_options.yml#trimming_options
            - "null"
    bqsr_known_sites:
        type: File[]
        secondaryFiles: [.tbi]
        doc: "One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis."
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
    qc_intervals:
        type: File
    picard_metric_accumulation_level:
        type: string
    roi_intervals:
       type: File
    vep_cache_dir:
        type:
            - string
            - Directory
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
    vep_pick:
        type:
            - "null"
            - type: enum
              symbols: ["pick", "flag_pick", "pick_allele", "per_gene", "pick_allele_gene", "flag_pick_allele", "flag_pick_allele_gene"]
    sample_name:
        type: string
    docm_vcf:
        type: File
        doc: "Common mutations in cancer that will be genotyped and passed through into the merged VCF if they have even low-level evidence of a mutation (by default, marked with filter DOCM_ONLY)"
    vep_custom_annotations:
        type: ../types/vep_custom_annotation.yml#vep_custom_annotation[]
        doc: "custom type, check types directory for input format"
    readcount_minimum_mapping_quality:
        type: int?
    readcount_minimum_base_quality:
        type: int?
    per_base_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    per_target_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    summary_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
outputs:
    cram:
        type: File
        outputSource: index_cram/indexed_cram
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
    tumor_snv_bam_readcount_tsv:
        type: File
        outputSource: detect_variants/tumor_snv_bam_readcount_tsv
    tumor_indel_bam_readcount_tsv:
        type: File
        outputSource: detect_variants/tumor_indel_bam_readcount_tsv
    per_base_coverage_metrics:
        type: File[]
        outputSource: alignment_and_qc/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File[]
        outputSource: alignment_and_qc/per_base_hs_metrics
    per_target_coverage_metrics:
        type: File[]
        outputSource: alignment_and_qc/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File[]
        outputSource: alignment_and_qc/per_target_hs_metrics
    summary_hs_metrics:
        type: File[]
        outputSource: alignment_and_qc/summary_hs_metrics
steps:
    alignment_and_qc:
        run: alignment_wgs.cwl
        in:
            reference: reference
            sequence: sequence
            trimming: trimming
            bqsr_known_sites: bqsr_known_sites
            omni_vcf: omni_vcf
            intervals: qc_intervals
            picard_metric_accumulation_level: picard_metric_accumulation_level
            minimum_mapping_quality: readcount_minimum_mapping_quality
            minimum_base_quality: readcount_minimum_base_quality
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
        out:
            [bam, mark_duplicates_metrics, insert_size_metrics, insert_size_histogram, alignment_summary_metrics, gc_bias_metrics, gc_bias_metrics_chart, gc_bias_metrics_summary, wgs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth, per_base_coverage_metrics, per_base_hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, summary_hs_metrics]
    detect_variants:
        run: tumor_only_detect_variants.cwl
        in:
            reference: reference
            bam: alignment_and_qc/bam
            roi_intervals: roi_intervals
            #varscan_strand_filter:
            #varscan_min_coverage:
            #varscan_min_var_freq:
            #varscan_p_value:
            #varscan_min_reads:
            #maximum_population_allele_frequency:
            vep_cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            vep_pick: vep_pick
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            #variants_to_table_fields:
            #variants_to_table_genotype_fields:
            #vep_to_table_fields:
            sample_name: sample_name
            docm_vcf: docm_vcf
            vep_custom_annotations: vep_custom_annotations
            readcount_minimum_mapping_quality: readcount_minimum_mapping_quality
            readcount_minimum_base_quality: readcount_minimum_base_quality
        out:
            [varscan_vcf, docm_gatk_vcf, annotated_vcf, final_vcf, final_tsv, vep_summary, tumor_snv_bam_readcount_tsv, tumor_indel_bam_readcount_tsv]
    bam_to_cram:
        run: ../tools/bam_to_cram.cwl
        in:
            bam: alignment_and_qc/bam
            reference: reference
        out:
            [cram]
    index_cram:
         run: ../tools/index_cram.cwl
         in:
            cram: bam_to_cram/cram
         out:
            [indexed_cram]
