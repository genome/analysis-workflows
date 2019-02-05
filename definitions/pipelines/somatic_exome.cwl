#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "exome alignment and somatic variant detection"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    reference: string
    tumor_bams:
        type: File[]
    tumor_readgroups:
        type: string[]
    tumor_cram_name:
        type: string?
        default: 'tumor.cram'
    normal_bams:
        type: File[]
    normal_readgroups:
        type: string[]
    normal_cram_name:
        type: string?
        default: 'normal.cram'
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
        type: string[]
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
    qc_minimum_mapping_quality:
        type: int?
        default: 0
    qc_minimum_base_quality:
        type: int?
        default: 0
    interval_list:
        type: File
    cosmic_vcf:
        type: File?
        secondaryFiles: [.tbi]
    panel_of_normals_vcf:
        type: File?
        secondaryFiles: [.tbi]
    strelka_cpu_reserved:
        type: int?
        default: 8
    mutect_scatter_count:
        type: int
    mutect_artifact_detection_mode:
        type: boolean
        default: false
    mutect_max_alt_allele_in_normal_fraction:
        type: float?
    mutect_max_alt_alleles_in_normal_count:
        type: int?
    varscan_strand_filter:
        type: int?
        default: 0
    varscan_min_coverage:
        type: int?
        default: 8
    varscan_min_var_freq:
        type: float?
        default: 0.05
    varscan_p_value:
        type: float?
        default: 0.99
    varscan_max_normal_freq:
        type: float?
    pindel_insert_size:
        type: int
        default: 400
    docm_vcf:
        type: File
        secondaryFiles: [.tbi]
    vep_cache_dir:
        type: string
    synonyms_file:
        type: File?
    annotate_coding_only:
        type: boolean?
    vep_pick:
        type:
            - "null"
            - type: enum
              symbols: ["pick", "flag_pick", "pick_allele", "per_gene", "pick_allele_gene", "flag_pick_allele", "flag_pick_allele_gene"]
    cle_vcf_filter:
        type: boolean
        default: false
    variants_to_table_fields:
        type: string[]
        default: [CHROM,POS,ID,REF,ALT,set,AC,AF]
    variants_to_table_genotype_fields:
        type: string[]
        default: [GT,AD]
    vep_to_table_fields:
        type: string[]
        default: [HGVSc,HGVSp]
    custom_gnomad_vcf:
        type: File?
        secondaryFiles: [.tbi]
    custom_clinvar_vcf:
        type: File?
        secondaryFiles: [.tbi]
outputs:
    tumor_cram:
        type: File
        outputSource: tumor_alignment_and_qc/cram
    tumor_mark_duplicates_metrics:
        type: File
        outputSource: tumor_alignment_and_qc/mark_duplicates_metrics
    tumor_insert_size_metrics:
        type: File
        outputSource: tumor_alignment_and_qc/insert_size_metrics
    tumor_alignment_summary_metrics:
        type: File
        outputSource: tumor_alignment_and_qc/alignment_summary_metrics
    tumor_hs_metrics:
        type: File
        outputSource: tumor_alignment_and_qc/hs_metrics
    tumor_per_target_coverage_metrics:
        type: File[]
        outputSource: tumor_alignment_and_qc/per_target_coverage_metrics
    tumor_per_target_hs_metrics:
        type: File[]
        outputSource: tumor_alignment_and_qc/per_target_hs_metrics
    tumor_per_base_coverage_metrics:
        type: File[]
        outputSource: tumor_alignment_and_qc/per_base_coverage_metrics
    tumor_per_base_hs_metrics:
        type: File[]
        outputSource: tumor_alignment_and_qc/per_base_hs_metrics
    tumor_summary_hs_metrics:
        type: File[]
        outputSource: tumor_alignment_and_qc/summary_hs_metrics
    tumor_flagstats:
        type: File
        outputSource: tumor_alignment_and_qc/flagstats
    tumor_verify_bam_id_metrics:
        type: File
        outputSource: tumor_alignment_and_qc/verify_bam_id_metrics
    tumor_verify_bam_id_depth:
        type: File
        outputSource: tumor_alignment_and_qc/verify_bam_id_depth
    normal_cram:
        type: File
        outputSource: normal_alignment_and_qc/cram
    normal_mark_duplicates_metrics:
        type: File
        outputSource: normal_alignment_and_qc/mark_duplicates_metrics
    normal_insert_size_metrics:
        type: File
        outputSource: normal_alignment_and_qc/insert_size_metrics
    normal_alignment_summary_metrics:
        type: File
        outputSource: normal_alignment_and_qc/alignment_summary_metrics
    normal_hs_metrics:
        type: File
        outputSource: normal_alignment_and_qc/hs_metrics
    normal_per_target_coverage_metrics:
        type: File[]
        outputSource: normal_alignment_and_qc/per_target_coverage_metrics
    normal_per_target_hs_metrics:
        type: File[]
        outputSource: normal_alignment_and_qc/per_target_hs_metrics
    normal_per_base_coverage_metrics:
        type: File[]
        outputSource: normal_alignment_and_qc/per_base_coverage_metrics
    normal_per_base_hs_metrics:
        type: File[]
        outputSource: normal_alignment_and_qc/per_base_hs_metrics
    normal_summary_hs_metrics:
        type: File[]
        outputSource: normal_alignment_and_qc/summary_hs_metrics
    normal_flagstats:
        type: File
        outputSource: normal_alignment_and_qc/flagstats
    normal_verify_bam_id_metrics:
        type: File
        outputSource: normal_alignment_and_qc/verify_bam_id_metrics
    normal_verify_bam_id_depth:
        type: File
        outputSource: normal_alignment_and_qc/verify_bam_id_depth
    mutect_unfiltered_vcf:
        type: File
        outputSource: detect_variants/mutect_unfiltered_vcf
        secondaryFiles: [.tbi]
    mutect_filtered_vcf:
        type: File
        outputSource: detect_variants/mutect_filtered_vcf
        secondaryFiles: [.tbi]
    strelka_unfiltered_vcf:
        type: File
        outputSource: detect_variants/strelka_unfiltered_vcf
        secondaryFiles: [.tbi]
    strelka_filtered_vcf:
        type: File
        outputSource: detect_variants/strelka_filtered_vcf
        secondaryFiles: [.tbi]
    varscan_unfiltered_vcf:
        type: File
        outputSource: detect_variants/varscan_unfiltered_vcf
        secondaryFiles: [.tbi]
    varscan_filtered_vcf:
        type: File
        outputSource: detect_variants/varscan_filtered_vcf
        secondaryFiles: [.tbi]
    pindel_unfiltered_vcf:
        type: File
        outputSource: detect_variants/pindel_unfiltered_vcf
        secondaryFiles: [.tbi]
    pindel_filtered_vcf:
        type: File
        outputSource: detect_variants/pindel_filtered_vcf
        secondaryFiles: [.tbi]
    docm_unfiltered_vcf:
        type: File
        outputSource: detect_variants/docm_unfiltered_vcf
    docm_filtered_vcf:
        type: File
        outputSource: detect_variants/docm_filtered_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: detect_variants/final_vcf
        secondaryFiles: [.tbi]
    final_filtered_vcf:
        type: File
        outputSource: detect_variants/final_filtered_vcf
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
    normal_snv_bam_readcount_tsv:
        type: File
        outputSource: detect_variants/normal_snv_bam_readcount_tsv
    normal_indel_bam_readcount_tsv:
        type: File
        outputSource: detect_variants/normal_indel_bam_readcount_tsv
    intervals_antitarget:
        type: File?
        outputSource: cnvkit/intervals_antitarget
    intervals_target:
        type: File?
        outputSource: cnvkit/intervals_target
    normal_antitarget_coverage:
        type: File
        outputSource: cnvkit/normal_antitarget_coverage
    normal_target_coverage:
        type: File
        outputSource: cnvkit/normal_target_coverage
    reference_coverage:
        type: File?
        outputSource: cnvkit/reference_coverage
    cn_diagram:
        type: File?
        outputSource: cnvkit/cn_diagram
    cn_scatter_plot:
        type: File?
        outputSource: cnvkit/cn_scatter_plot
    tumor_antitarget_coverage:
        type: File
        outputSource: cnvkit/tumor_antitarget_coverage
    tumor_target_coverage:
        type: File
        outputSource: cnvkit/tumor_target_coverage
    tumor_bin_level_ratios:
        type: File
        outputSource: cnvkit/tumor_bin_level_ratios
    tumor_segmented_ratios:
        type: File
        outputSource: cnvkit/tumor_segmented_ratios
steps:
    tumor_alignment_and_qc:
        run: exome_alignment.cwl
        in:
            reference: reference
            bams: tumor_bams
            readgroups: tumor_readgroups
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
            minimum_mapping_quality: qc_minimum_mapping_quality
            minimum_base_quality: qc_minimum_base_quality
            final_name: tumor_cram_name
        out:
            [cram, mark_duplicates_metrics, insert_size_metrics, alignment_summary_metrics, hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, per_base_coverage_metrics, per_base_hs_metrics, summary_hs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
    normal_alignment_and_qc:
        run: exome_alignment.cwl
        in:
            reference: reference
            bams: normal_bams
            readgroups: normal_readgroups
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
            minimum_mapping_quality: qc_minimum_mapping_quality
            minimum_base_quality: qc_minimum_base_quality
            final_name: normal_cram_name
        out:
            [cram, mark_duplicates_metrics, insert_size_metrics, alignment_summary_metrics, hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, per_base_coverage_metrics, per_base_hs_metrics, summary_hs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
    detect_variants:
        run: detect_variants.cwl
        in:
            reference: reference
            tumor_cram: tumor_alignment_and_qc/cram
            normal_cram: normal_alignment_and_qc/cram
            interval_list: interval_list
            dbsnp_vcf: dbsnp_vcf
            cosmic_vcf: cosmic_vcf
            panel_of_normals_vcf: panel_of_normals_vcf
            strelka_exome_mode:
                default: true
            strelka_cpu_reserved: strelka_cpu_reserved
            mutect_scatter_count: mutect_scatter_count
            mutect_artifact_detection_mode: mutect_artifact_detection_mode
            mutect_max_alt_allele_in_normal_fraction: mutect_max_alt_allele_in_normal_fraction
            mutect_max_alt_alleles_in_normal_count: mutect_max_alt_alleles_in_normal_count
            varscan_strand_filter: varscan_strand_filter
            varscan_min_coverage: varscan_min_coverage
            varscan_min_var_freq: varscan_min_var_freq
            varscan_p_value: varscan_p_value
            varscan_max_normal_freq: varscan_max_normal_freq
            pindel_insert_size: pindel_insert_size
            docm_vcf: docm_vcf
            vep_cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            vep_pick: vep_pick
            cle_vcf_filter: cle_vcf_filter
            variants_to_table_fields: variants_to_table_fields
            variants_to_table_genotype_fields: variants_to_table_genotype_fields
            vep_to_table_fields: vep_to_table_fields
            custom_gnomad_vcf: custom_gnomad_vcf
            custom_clinvar_vcf: custom_clinvar_vcf
        out:
            [mutect_unfiltered_vcf, mutect_filtered_vcf, strelka_unfiltered_vcf, strelka_filtered_vcf, varscan_unfiltered_vcf, varscan_filtered_vcf, pindel_unfiltered_vcf, pindel_filtered_vcf, docm_unfiltered_vcf, docm_filtered_vcf, final_vcf, final_filtered_vcf, final_tsv, vep_summary, tumor_snv_bam_readcount_tsv, tumor_indel_bam_readcount_tsv, normal_snv_bam_readcount_tsv, normal_indel_bam_readcount_tsv]
    cnvkit:
        run: ../subworkflows/cram_to_cnvkit.cwl
        in: 
            normal_cram: tumor_alignment_and_qc/cram
            tumor_cram: normal_alignment_and_qc/cram
            reference: reference
            bait_intervals: bait_intervals
        out:
            [intervals_antitarget, intervals_target, normal_antitarget_coverage, normal_target_coverage, reference_coverage, cn_diagram, cn_scatter_plot, tumor_antitarget_coverage, tumor_target_coverage, tumor_bin_level_ratios, tumor_segmented_ratios]
