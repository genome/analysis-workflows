#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Detect Variants workflow"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type: string
    tumor_bam:
        type: File
        secondaryFiles: [.bai,^.bai]
    normal_bam:
        type: File
        secondaryFiles: [.bai,^.bai]
    interval_list:
        type: File
    strelka_exome_mode:
        type: boolean
    strelka_cpu_reserved:
        type: int?
        default: 8
    readcount_minimum_base_quality:
        type: int?
    readcount_minimum_mapping_quality:
        type: int?
    mutect_scatter_count:
        type: int?
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
    varscan_max_normal_freq:
        type: float?
    pindel_insert_size:
        type: int
        default: 400
    docm_vcf:
        type: File
        secondaryFiles: [.tbi]
    filter_docm_variants:
        type: boolean?
        default: true
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
    vep_pick:
        type:
            - "null"
            - type: enum
              symbols: ["pick", "flag_pick", "pick_allele", "per_gene", "pick_allele_gene", "flag_pick_allele", "flag_pick_allele_gene"]
    vep_plugins:
        type: string[]?
        default: [Downstream, Wildtype]
    filter_gnomADe_maximum_population_allele_frequency:
        type: float?
        default: 0.001
    filter_mapq0_threshold:
        type: float?
        default: 0.15
    filter_minimum_depth:
        type: int?
        default: 1
    cle_vcf_filter:
        type: boolean?
        default: false
    filter_somatic_llr_threshold:
        type: float?
        default: 5
    variants_to_table_fields:
        type: string[]?
        default: [CHROM,POS,ID,REF,ALT,set,AC,AF]
    variants_to_table_genotype_fields:
        type: string[]?
        default: [GT,AD]
    vep_to_table_fields:
        type: string[]?
        default: [HGVSc,HGVSp]
    custom_gnomad_vcf:
        type: File?
        secondaryFiles: [.tbi]
    custom_clinvar_vcf:
        type: File?
        secondaryFiles: [.tbi]
outputs:
    mutect_unfiltered_vcf:
        type: File
        outputSource: mutect/unfiltered_vcf
        secondaryFiles: [.tbi]
    mutect_filtered_vcf:
        type: File
        outputSource: mutect/filtered_vcf
        secondaryFiles: [.tbi]
    strelka_unfiltered_vcf:
        type: File
        outputSource: strelka/unfiltered_vcf
        secondaryFiles: [.tbi]
    strelka_filtered_vcf:
        type: File
        outputSource: strelka/filtered_vcf
        secondaryFiles: [.tbi]
    varscan_unfiltered_vcf:
        type: File
        outputSource: varscan/unfiltered_vcf
        secondaryFiles: [.tbi]
    varscan_filtered_vcf:
        type: File
        outputSource: varscan/filtered_vcf
        secondaryFiles: [.tbi]
    pindel_unfiltered_vcf:
        type: File
        outputSource: pindel/unfiltered_vcf
        secondaryFiles: [.tbi]
    pindel_filtered_vcf:
        type: File
        outputSource: pindel/filtered_vcf
        secondaryFiles: [.tbi]
    docm_filtered_vcf:
        type: File
        outputSource: docm/docm_variants_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: rcnt_and_filter_variants/final_vcf
        secondaryFiles: [.tbi]
    final_filtered_vcf:
        type: File
        outputSource: rcnt_and_filter_variants/final_filtered_vcf
        secondaryFiles: [.tbi]
    final_tsv:
        type: File
        outputSource: rcnt_and_filter_variants/final_tsv
    vep_summary:
        type: File
        outputSource: annotate_variants/vep_summary
    tumor_snv_bam_readcount_tsv:
        type: File
        outputSource: rcnt_and_filter_variants/tumor_snv_bam_readcount_tsv
    tumor_indel_bam_readcount_tsv:
        type: File
        outputSource: rcnt_and_filter_variants/tumor_indel_bam_readcount_tsv
    normal_snv_bam_readcount_tsv:
        type: File
        outputSource: rcnt_and_filter_variants/normal_snv_bam_readcount_tsv
    normal_indel_bam_readcount_tsv:
        type: File
        outputSource: rcnt_and_filter_variants/normal_indel_bam_readcount_tsv
steps:
    mutect:
        run: ../subworkflows/mutect.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: interval_list
            scatter_count: mutect_scatter_count
        out:
            [unfiltered_vcf, filtered_vcf]
    strelka:
        run: ../subworkflows/strelka_and_post_processing.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: interval_list
            exome_mode: strelka_exome_mode
            cpu_reserved: strelka_cpu_reserved
        out:
            [unfiltered_vcf, filtered_vcf]
    varscan:
        run: ../subworkflows/varscan_pre_and_post_processing.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: interval_list
            strand_filter: varscan_strand_filter
            min_coverage: varscan_min_coverage
            min_var_freq: varscan_min_var_freq
            p_value: varscan_p_value
            max_normal_freq: varscan_max_normal_freq
        out:
            [unfiltered_vcf, filtered_vcf]
    pindel:
        run: ../subworkflows/pindel.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: interval_list
            insert_size: pindel_insert_size
        out:
            [unfiltered_vcf, filtered_vcf]
    docm:
        run: ../subworkflows/docm_cle.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            docm_vcf: docm_vcf
            interval_list: interval_list
            filter_docm_variants: filter_docm_variants
        out:
            [docm_variants_vcf]
    combine:
        run: ../tools/combine_variants.cwl
        in:
            reference: reference
            mutect_vcf: mutect/filtered_vcf
            strelka_vcf: strelka/filtered_vcf
            varscan_vcf: varscan/filtered_vcf
            pindel_vcf: pindel/filtered_vcf
        out:
            [combined_vcf]
    add_docm_variants:
        run: ../tools/docm_add_variants.cwl
        in:
            reference: reference
            docm_vcf: docm/docm_variants_vcf
            callers_vcf: combine/combined_vcf
        out:
            [merged_vcf]
    decompose:
        run: ../tools/vt_decompose.cwl
        in:
            vcf: add_docm_variants/merged_vcf
        out:
            [decomposed_vcf]
    decompose_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: decompose/decomposed_vcf
        out:
            [indexed_vcf]
    annotate_variants:
        run: ../tools/vep.cwl
        in:
            vcf: decompose_index/indexed_vcf
            cache_dir: vep_cache_dir
            ensembl_assembly: vep_ensembl_assembly
            ensembl_version: vep_ensembl_version
            ensembl_species: vep_ensembl_species
            synonyms_file: synonyms_file
            coding_only: annotate_coding_only
            reference: reference
            custom_gnomad_vcf: custom_gnomad_vcf
            pick: vep_pick
            custom_clinvar_vcf: custom_clinvar_vcf
            plugins: vep_plugins
        out:
            [annotated_vcf, vep_summary]
    rcnt_and_filter_variants:
        run: ../subworkflows/rcnt_and_filter_vcf.cwl
        in:
            annotated_vcf: annotate_variants/annotated_vcf
            reference: reference
            tumor_cram: tumor_bam
            normal_cram: normal_bam
            readcount_minimum_base_quality: readcount_minimum_base_quality
            readcount_minimum_mapping_quality: readcount_minimum_mapping_quality
            filter_gnomADe_maximum_population_allele_frequency: filter_gnomADe_maximum_population_allele_frequency
            filter_mapq0_threshold: filter_mapq0_threshold
            filter_minimum_depth: filter_minimum_depth
            cle_vcf_filter: cle_vcf_filter
            filter_somatic_llr_threshold: filter_somatic_llr_threshold
            variants_to_table_fields: variants_to_table_fields
            variants_to_table_genotype_fields: variants_to_table_genotype_fields
            vep_to_table_fields: variants_to_table_fields
        out: [ final_vcf, final_filtered_vcf, final_tsv, tumor_snv_bam_readcount_tsv, tumor_indel_bam_readcount_tsv, normal_snv_bam_readcount_tsv, normal_indel_bam_readcount_tsv]
