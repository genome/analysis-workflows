#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Sub workflow for rcnt and filtering of variants"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    annotated_vcf:
        type: File
    reference:
        type: string
    tumor_cram:
        type: File
        secondaryFiles: [.crai,^.crai]
    normal_cram:
        type: File
        secondaryFiles: [.crai,^.crai]
    readcount_minimum_base_quality:
        type: int?
    readcount_minimum_mapping_quality:
        type: int?
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
outputs:
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
    tumor_cramToBam:
        run: ../tools/samtools_view_convert2bam.cwl
        in:
            file: tumor_cram
        out: [ bam_file ]
    tumor_indexBam:
        run: ../tools/index_bam.cwl
        in:
            bam: tumor_cramToBam/bam_file
        out: [indexed_bam]
    normal_cramToBam:
        run: ../tools/samtools_view_convert2bam.cwl
        in:
            file: normal_cram
        out: [ bam_file ]
    normal_indexBam:
        run: ../tools/index_bam.cwl
        in:
            bam: normal_cramToBam/bam_file
        out: [ indexed_bam ]
    tumor_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: annotated_vcf
            sample:
                default: 'TUMOR'
            reference_fasta: reference
            bam: tumor_indexBam/indexed_bam
            min_base_quality: readcount_minimum_base_quality
            min_mapping_quality: readcount_minimum_mapping_quality
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    normal_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: annotated_vcf
            sample:
                default: 'NORMAL'
            reference_fasta: reference
            bam: normal_indexBam/indexed_bam
            min_base_quality: readcount_minimum_base_quality
            min_mapping_quality: readcount_minimum_mapping_quality
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    add_tumor_bam_readcount_to_vcf:
        run: ../subworkflows/vcf_readcount_annotator.cwl
        in:
            vcf: annotated_vcf
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
            filter_minimum_depth: filter_minimum_depth
            sample_names:
                default: 'NORMAL,TUMOR'
            tumor_bam: tumor_indexBam/indexed_bam
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
