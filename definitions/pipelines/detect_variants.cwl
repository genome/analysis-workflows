#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Detect Variants workflow"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type: string
    tumor_cram:
        type: File
        secondaryFiles: [.crai,^.crai]
    normal_cram:
        type: File
        secondaryFiles: [.crai,^.crai]
    interval_list:
        type: File
    dbsnp_vcf:
        type: File?
        secondaryFiles: [.tbi]
    cosmic_vcf:
        type: File?
        secondaryFiles: [.tbi]
    panel_of_normals_vcf:
        type: File?
        secondaryFiles: [.tbi]
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
    mutect_artifact_detection_mode:
        type: boolean?
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
    filter_gnomADe_maximum_population_allele_frequency:
        type: float?
        default: 0.001
    filter_mapq0_threshold:
        type: float?
        default: 0.15
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
    docm_unfiltered_vcf:
        type: File
        outputSource: docm/unfiltered_vcf
    docm_filtered_vcf:
        type: File
        outputSource: docm/filtered_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [.tbi]
    final_filtered_vcf:
        type: File
        outputSource: annotated_filter_index/indexed_vcf
        secondaryFiles: [.tbi]
    final_tsv:
        type: File
        outputSource: add_vep_fields_to_table/annotated_variants_tsv
    vep_summary:
        type: File
        outputSource: annotate_variants/vep_summary
    tumor_snv_bam_readcount_tsv:
        type: File
        outputSource: tumor_bam_readcount/snv_bam_readcount_tsv
    tumor_indel_bam_readcount_tsv:
        type: File
        outputSource: tumor_bam_readcount/indel_bam_readcount_tsv
    normal_snv_bam_readcount_tsv:
        type: File
        outputSource: normal_bam_readcount/snv_bam_readcount_tsv
    normal_indel_bam_readcount_tsv:
        type: File
        outputSource: normal_bam_readcount/indel_bam_readcount_tsv
steps:
    mutect:
        run: ../subworkflows/mutect.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            interval_list: interval_list
            dbsnp_vcf: dbsnp_vcf
            cosmic_vcf: cosmic_vcf
            max_alt_allele_in_normal_fraction: mutect_max_alt_allele_in_normal_fraction
            max_alt_alleles_in_normal_count: mutect_max_alt_alleles_in_normal_count
            scatter_count: mutect_scatter_count
            artifact_detection_mode: mutect_artifact_detection_mode
            panel_of_normals_vcf: panel_of_normals_vcf
        out:
            [unfiltered_vcf, filtered_vcf]
    strelka:
        run: ../subworkflows/strelka_and_post_processing.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            interval_list: interval_list
            exome_mode: strelka_exome_mode
            cpu_reserved: strelka_cpu_reserved
        out:
            [unfiltered_vcf, filtered_vcf]
    varscan:
        run: ../subworkflows/varscan_pre_and_post_processing.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
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
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            interval_list: interval_list
            insert_size: pindel_insert_size
        out:
            [unfiltered_vcf, filtered_vcf]
    docm:
        run: ../subworkflows/docm_cle.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            docm_vcf: docm_vcf
            interval_list: interval_list
        out:
            [unfiltered_vcf, filtered_vcf]
    combine:
        run: ../tools/combine_variants.cwl
        in:
            reference: reference
            mutect_vcf: mutect/filtered_vcf
            strelka_vcf: strelka/filtered_vcf
            varscan_vcf: varscan/filtered_vcf
            pindel_vcf: pindel/filtered_vcf
            docm_vcf: docm/filtered_vcf
        out:
            [combined_vcf]
    decompose:
        run: ../tools/vt_decompose.cwl
        in:
            vcf: combine/combined_vcf
        out:
            [decomposed_vcf]
    annotate_variants:
        run: ../tools/vep.cwl
        in:
            vcf: decompose/decomposed_vcf
            cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            coding_only: annotate_coding_only
            reference: reference
            custom_gnomad_vcf: custom_gnomad_vcf
            pick: vep_pick
            custom_clivnar_vcf: custom_clinvar_vcf
        out:
            [annotated_vcf, vep_summary]
    tumor_cram_to_bam:
        run: ../subworkflows/cram_to_bam_and_index.cwl
        in:
            cram: tumor_cram
            reference: reference
        out:
            [bam]
    normal_cram_to_bam:
        run: ../subworkflows/cram_to_bam_and_index.cwl
        in:
            cram: normal_cram
            reference: reference
        out:
            [bam]
    tumor_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: annotate_variants/annotated_vcf
            sample:
                default: 'TUMOR'
            reference_fasta: reference
            bam: tumor_cram_to_bam/bam
            min_base_quality: readcount_minimum_base_quality
            min_mapping_quality: readcount_minimum_mapping_quality
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    normal_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: annotate_variants/annotated_vcf
            sample:
                default: 'NORMAL'
            reference_fasta: reference
            bam: normal_cram_to_bam/bam
            min_base_quality: readcount_minimum_base_quality
            min_mapping_quality: readcount_minimum_mapping_quality
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    add_tumor_bam_readcount_to_vcf:
        run: ../subworkflows/vcf_readcount_annotator.cwl
        in:
            vcf: annotate_variants/annotated_vcf
            snv_bam_readcount_tsv: tumor_bam_readcount/snv_bam_readcount_tsv
            indel_bam_readcount_tsv: tumor_bam_readcount/indel_bam_readcount_tsv
            data_type:
                default: 'DNA'
            sample_name:
                default: 'TUMOR'
        out:
            [annotated_bam_readcount_vcf]
    add_normal_bam_readcount_to_vcf:
        run: ../subworkflows/vcf_readcount_annotator.cwl
        in:
            vcf: add_tumor_bam_readcount_to_vcf/annotated_bam_readcount_vcf
            snv_bam_readcount_tsv: normal_bam_readcount/snv_bam_readcount_tsv
            indel_bam_readcount_tsv: normal_bam_readcount/indel_bam_readcount_tsv
            data_type:
                default: 'DNA'
            sample_name:
                default: 'NORMAL'
        out:
            [annotated_bam_readcount_vcf]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: add_normal_bam_readcount_to_vcf/annotated_bam_readcount_vcf
        out:
            [indexed_vcf]
    filter_vcf:
        run: ../subworkflows/filter_vcf.cwl
        in: 
            vcf: index/indexed_vcf
            filter_gnomADe_maximum_population_allele_frequency: filter_gnomADe_maximum_population_allele_frequency
            filter_mapq0_threshold: filter_mapq0_threshold
            filter_somatic_llr_threshold: filter_somatic_llr_threshold
            tumor_bam: tumor_cram_to_bam/bam
            do_cle_vcf_filter: cle_vcf_filter
            reference: reference
        out: 
            [filtered_vcf]
    annotated_filter_bgzip:
        run: ../tools/bgzip.cwl
        in:
            file: filter_vcf/filtered_vcf
        out:
            [bgzipped_file]
    annotated_filter_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: annotated_filter_bgzip/bgzipped_file
        out:
            [indexed_vcf]
    variants_to_table:
        run: ../tools/variants_to_table.cwl
        in:
            reference: reference
            vcf: annotated_filter_index/indexed_vcf
            fields: variants_to_table_fields
            genotype_fields: variants_to_table_genotype_fields
        out:
            [variants_tsv]
    add_vep_fields_to_table:
        run: ../tools/add_vep_fields_to_table.cwl
        in:
            vcf: annotated_filter_index/indexed_vcf
            vep_fields: vep_to_table_fields
            tsv: variants_to_table/variants_tsv
        out: [annotated_variants_tsv]
