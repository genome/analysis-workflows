#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "exome alignment and tumor-only variant detection"
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
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa, .index]
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
    varscan_strand_filter:
        type: int?
        default: 0
    varscan_min_coverage:
        type: int?
        default: 8
    varscan_min_var_freq:
        type: float?
        default: 0.1
    varscan_p_value:
        type: float?
        default: 0.99
    varscan_min_reads:
        type: int?
        default: 2
    maximum_population_allele_frequency:
        type: float?
        default: 0.001
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
        default: true
    vep_pick:
        type:
            - "null"
            - type: enum
              symbols: ["pick", "flag_pick", "pick_allele", "per_gene", "pick_allele_gene", "flag_pick_allele", "flag_pick_allele_gene"]
    variants_to_table_fields:
        type: string[]?
        default: [CHROM,POS,REF,ALT,set]
    variants_to_table_genotype_fields:
        type: string[]?
        default: [GT,AD,AF,DP]
    vep_to_table_fields:
        type: string[]?
        default: [Consequence,SYMBOL,Feature_type,Feature,HGVSc,HGVSp,cDNA_position,CDS_position,Protein_position,Amino_acids,Codons,HGNC_ID,Existing_variation,gnomADe_AF,CLIN_SIG,SOMATIC,PHENO]
    sample_name:
        type: string
    docm_vcf:
        type: File
        secondaryFiles: [.tbi]
    vep_custom_annotations:
        type: ../types/vep_custom_annotation.yml#vep_custom_annotation[]
        doc: "custom type, check types directory for input format"
    qc_minimum_mapping_quality:
        type: int?
    qc_minimum_base_quality:
        type: int?
    readcount_minimum_mapping_quality:
        type: int?
    readcount_minimum_base_quality:
        type: int?
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
    hs_metrics:
        type: File
        outputSource: alignment_and_qc/hs_metrics
    per_target_coverage_metrics:
        type: File[]
        outputSource: alignment_and_qc/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File[]
        outputSource: alignment_and_qc/per_target_hs_metrics
    per_base_coverage_metrics:
        type: File[]
        outputSource: alignment_and_qc/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File[]
        outputSource: alignment_and_qc/per_base_hs_metrics
    summary_hs_metrics:
        type: File[]
        outputSource: alignment_and_qc/summary_hs_metrics
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
steps:
    alignment_and_qc:
        run: alignment_exome.cwl
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
            qc_minimum_mapping_quality: qc_minimum_mapping_quality
            qc_minimum_base_quality: qc_minimum_base_quality
        out:
            [bam, mark_duplicates_metrics, insert_size_metrics, insert_size_histogram, alignment_summary_metrics, hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, per_base_coverage_metrics, per_base_hs_metrics, summary_hs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
    detect_variants:
        run: tumor_only_detect_variants.cwl
        in:
            reference: reference
            bam: alignment_and_qc/bam
            interval_list: target_intervals
            varscan_strand_filter: varscan_strand_filter
            varscan_min_coverage: varscan_min_coverage
            varscan_min_var_freq: varscan_min_var_freq
            varscan_p_value: varscan_p_value
            varscan_min_reads: varscan_min_reads
            maximum_population_allele_frequency: maximum_population_allele_frequency
            vep_cache_dir: vep_cache_dir
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            synonyms_file: synonyms_file
            vep_pick: vep_pick
            variants_to_table_fields: variants_to_table_fields
            variants_to_table_genotype_fields: variants_to_table_genotype_fields
            vep_to_table_fields: vep_to_table_fields
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
