#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Tumor-Only Detect Variants workflow"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement
inputs:
    reference:
        type: string
    cram:
        type: File
        secondaryFiles: [^.crai,.crai]
    interval_list:
        type: File
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
    synonyms_file:
        type: File?
    coding_only:
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
    custom_gnomad_vcf:
        type: File?
        secondaryFiles: [.tbi]
    readcount_minimum_mapping_quality:
        type: int?
    readcount_minimum_base_quality:
        type: int?
outputs:
    varscan_vcf:
        type: File
        outputSource: varscan/unfiltered_vcf
        secondaryFiles: [.tbi]
    docm_gatk_vcf:
        type: File
        outputSource: docm/unfiltered_vcf
    annotated_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: index_filtered/indexed_vcf
        secondaryFiles: [.tbi]
    final_tsv:
        type: File
        outputSource: add_vep_fields_to_table/annotated_variants_tsv
    vep_summary:
        type: File
        outputSource: annotate_variants/vep_summary
    tumor_snv_bam_readcount_tsv:
        type: File
        outputSource: bam_readcount/snv_bam_readcount_tsv
    tumor_indel_bam_readcount_tsv:
        type: File
        outputSource: bam_readcount/indel_bam_readcount_tsv
steps:
    varscan:
        run: ../subworkflows/varscan_germline.cwl
        in:
            reference: reference
            cram: cram
            interval_list: interval_list
            strand_filter: varscan_strand_filter
            min_coverage: varscan_min_coverage
            min_var_freq: varscan_min_var_freq
            min_reads: varscan_min_reads
            p_value: varscan_p_value
            sample_name: sample_name
        out:
            [unfiltered_vcf, filtered_vcf]
    docm:
        run: ../subworkflows/docm_germline.cwl
        in:
            reference: reference
            cram: cram
            interval_list: interval_list
            docm_vcf: docm_vcf
        out:
            [unfiltered_vcf, filtered_vcf]
    combine_variants:
        run: ../tools/germline_combine_variants.cwl
        in:
            reference: reference
            varscan_vcf: varscan/filtered_vcf
            docm_vcf: docm/filtered_vcf
        out:
            [combined_vcf]
    decompose:
        run: ../tools/vt_decompose.cwl
        in:
            vcf: combine_variants/combined_vcf
        out:
            [decomposed_vcf]
    annotate_variants:
        run: ../tools/vep.cwl
        in:
            vcf: decompose/decomposed_vcf
            cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            coding_only: coding_only
            reference: reference
            custom_gnomad_vcf: custom_gnomad_vcf
            pick: vep_pick
        out:
            [annotated_vcf, vep_summary]
    cram_to_bam:
        run: ../subworkflows/cram_to_bam_and_index.cwl
        in:
            cram: cram
            reference: reference
        out:
            [bam]
    bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: annotate_variants/annotated_vcf
            sample: sample_name
            reference_fasta: reference
            bam: cram_to_bam/bam
            min_mapping_quality: readcount_minimum_mapping_quality
            min_base_quality: readcount_minimum_base_quality
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    add_bam_readcount_to_vcf:
        run: ../subworkflows/vcf_readcount_annotator.cwl
        in:
            vcf: annotate_variants/annotated_vcf
            snv_bam_readcount_tsv: bam_readcount/snv_bam_readcount_tsv
            indel_bam_readcount_tsv: bam_readcount/indel_bam_readcount_tsv
            data_type:
                default: 'DNA'
            sample_name: sample_name
        out:
            [annotated_bam_readcount_vcf]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: add_bam_readcount_to_vcf/annotated_bam_readcount_vcf
        out:
            [indexed_vcf]
    hard_filter:
        run: ../tools/select_variants.cwl
        in:
            reference: reference
            vcf: index/indexed_vcf
            interval_list: interval_list
            exclude_filtered:
                default: true
        out:
            [filtered_vcf]
    af_filter:
        run: ../tools/filter_vcf_gnomADe_allele_freq.cwl
        in:
            vcf: hard_filter/filtered_vcf
            maximum_population_allele_frequency: maximum_population_allele_frequency
        out:
            [filtered_vcf]
    coding_variant_filter:
        run: ../tools/filter_vcf_coding_variant.cwl
        in:
            vcf: af_filter/filtered_vcf
        out:
            [filtered_vcf]
    bgzip_filtered:
        run: ../tools/bgzip.cwl
        in:
            file: coding_variant_filter/filtered_vcf
        out:
            [bgzipped_file]
    index_filtered:
        run: ../tools/index_vcf.cwl
        in:
            vcf: bgzip_filtered/bgzipped_file
        out:
            [indexed_vcf]
    variants_to_table:
        run: ../tools/variants_to_table.cwl
        in:
            reference: reference
            vcf: index_filtered/indexed_vcf
            fields: variants_to_table_fields
            genotype_fields: variants_to_table_genotype_fields
        out:
            [variants_tsv]
    add_vep_fields_to_table:
        run: ../tools/add_vep_fields_to_table.cwl
        in:
            vcf: index_filtered/indexed_vcf
            vep_fields: vep_to_table_fields
            tsv: variants_to_table/variants_tsv
        out: [annotated_variants_tsv]
