#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "joint variant detection(snps,svs)"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/vep_custom_annotation.yml
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement
    - class: ScatterFeatureRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    bams:
        type: File[]
        secondaryFiles: [^.bai]
    sample_names:
        type: string[]
    cohort_name:
        type: string
    gvcf_gq_bands:
        type: string[]
    intervals:
        type:
            type: array
            items:
                type: array
                items: string
    contamination_fraction:
        type: string[]
    ploidy:
        type: int?
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
    vep_plugins:
        type: string[]
        default: [Frameshift, Wildtype]
    synonyms_file:
        type: File?
    annotate_coding_only:
        type: boolean?
    vep_custom_annotations:
        type: ../types/vep_custom_annotation.yml#vep_custom_annotation[]
        doc: "custom type, check types directory for input format"
    limit_variant_intervals:
        type: File
    snp_to_table_fields:
        type: string[]
        default: ['CHROM','POS','ID','REF','ALT']
    snp_to_table_genotype_fields:
        type: string[]
    vep_to_table_fields:
        type: string[]
    snp_final_tsv_prefix:
        type: string?
        default: 'variants'
    snp_gnomad_max_pop_af:
        type: float
        default: 0.05
    gatk_min_conf_call:
        type: float?


    sv_exclude_regions:
        type: File?
    manta_call_regions:
        type: File?
    manta_output_contigs:
        type: boolean?
    cnvnator_bin_size:
        type: int?
    cnvkit_method:
        type:
          - "null"
          - type: enum
            symbols: ["hybrid", "amplicon", "wgs"]
    cnvkit_reference_cnn:
        type: File?
    cnvkit_segment_filter:
        type:
          - "null"
          - type: enum
            symbols: ["ampdel", "ci", "cn", "sem"]
    sv_filter_del_depth:
        type: double?
    sv_filter_dup_depth:
        type: double?
    sv_filter_paired_count:
        type: int?
    sv_filter_split_count:
        type: int?
    sv_filter_alt_abundance_percentage:
        type: double?
    sv_filter_depth_caller_min_size:
        type: int?
    survivor_estimate_sv_distance:
        type: boolean
    survivor_max_distance_to_merge:
        type: int
    survivor_minimum_sv_calls:
        type: int
    survivor_minimum_sv_size:
        type: int
    survivor_same_strand:
        type: boolean
    survivor_same_type:
        type: boolean
    sv_filter_blocklist_bedpe:
        type: File?
    annotsv_filter_pop_af:
        type: double?
    annotsv_filter_no_CDS:
        type: boolean?
    annotsv_annotations:
        type:
            - string
            - Directory
outputs:
    snps_staged:
        type: Directory
        outputSource: detect_snps/all_staged
    svs_staged:
        type: Directory
        outputSource: detect_svs/all_staged
steps:
    detect_snps:
        run: joint_detect_snps.cwl
        in:
            reference: reference
            bams: bams
            sample_names: sample_names
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            contamination_fraction: contamination_fraction
            ploidy: ploidy
            vep_cache_dir: vep_cache_dir
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            vep_plugins: vep_plugins
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            vep_custom_annotations: vep_custom_annotations
            limit_variant_intervals: limit_variant_intervals
            variants_to_table_fields: snp_to_table_fields
            variants_to_table_genotype_fields: snp_to_table_genotype_fields
            vep_to_table_fields: vep_to_table_fields
            final_tsv_prefix: snp_final_tsv_prefix
            gnomad_max_pop_af: snp_gnomad_max_pop_af
            min_conf_call: gatk_min_conf_call
        out:
            [raw_vcf, all_staged]
    detect_svs:
        run: joint_detect_svs.cwl
        in:
            reference: reference
            bams: bams
            sample_names: sample_names
            cohort_name: cohort_name
            genome_build: vep_ensembl_assembly
            exclude_regions: sv_exclude_regions
            manta_call_regions: manta_call_regions
            manta_output_contigs: manta_output_contigs
            cnvnator_bin_size: cnvnator_bin_size
            cnvkit_method: cnvkit_method
            cnvkit_reference_cnn: cnvkit_reference_cnn
            cnvkit_segment_filter: cnvkit_segment_filter
            filter_del_depth: sv_filter_del_depth
            filter_dup_depth: sv_filter_dup_depth
            filter_paired_count: sv_filter_paired_count
            filter_split_count: sv_filter_split_count
            filter_alt_abundance_percentage: sv_filter_alt_abundance_percentage
            filter_depth_caller_min_size: sv_filter_depth_caller_min_size
            survivor_estimate_sv_distance: survivor_estimate_sv_distance
            survivor_max_distance_to_merge: survivor_max_distance_to_merge
            survivor_minimum_sv_calls: survivor_minimum_sv_calls
            survivor_minimum_sv_size: survivor_minimum_sv_size
            survivor_same_strand: survivor_same_strand
            survivor_same_type: survivor_same_type
            snps_vcf: detect_snps/raw_vcf
            filter_blocklist_bedpe: sv_filter_blocklist_bedpe
            annotsv_filter_pop_af: annotsv_filter_pop_af
            annotsv_annotations: annotsv_annotations
        out:
            [all_staged]
