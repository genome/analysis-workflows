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
    vep_cache_dir:
        type: string?
    synonyms_file:
        type: File?
    coding_only:
        type: boolean?
    hgvs_annotation:
        type: boolean?
    variants_to_table_fields:
        type: string[]?
        default: [CHROM,POS,ID,REF,ALT,set,AC,AF]
    variants_to_table_genotype_fields:
        type: string[]?
        default: [GT,AD]
    vep_to_table_fields:
        type: string[]?
        default: [HGVSc,HGVSp]
    sample_name:
        type: string
outputs:
    varscan_vcf:
        type: File
        outputSource: varscan/variants
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [.tbi]
    final_tsv:
        type: File
        outputSource: add_vep_fields_to_table/annotated_variants_tsv
    vep_summary:
        type: File
        outputSource: annotate_variants/vep_summary
    tumor_bam_readcount_tsv:
        type: File
        outputSource: bam_readcount/bam_readcount_tsv
steps:
    varscan:
        run: ../varscan/germline_workflow.cwl
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
            [varscan_vcf, variants]
    annotate_variants:
        run: vep.cwl
        in:
            vcf: varscan/varscan_vcf
            cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            coding_only: coding_only
            hgvs: hgvs_annotation
            reference: reference
        out:
            [annotated_vcf, vep_summary]
    cram_to_bam:
        run: ../cram_to_bam/workflow.cwl
        in:
            cram: cram
            reference: reference
        out:
            [bam]
    bam_readcount:
        run: ../pvacseq/bam_readcount.cwl
        in:
            vcf: varscan/varscan_vcf
            sample: sample_name
            reference_fasta: reference
            bam: cram_to_bam/bam
        out:
            [bam_readcount_tsv]
    add_bam_readcount_to_vcf:
        run: add_bam_readcount_to_vcf.cwl
        in:
            - id: vcf
              source: annotate_variants/annotated_vcf
            - id: bam_readcount_tsvs
              source: bam_readcount/bam_readcount_tsv
              valueFrom: ${ return [ self ]; }
            - id: sample_names
              source: sample_name
              valueFrom: ${ return [ self ]; }
        out:
            [annotated_bam_readcount_vcf]
    bgzip:
        run: bgzip.cwl
        in:
            file: add_bam_readcount_to_vcf/annotated_bam_readcount_vcf
        out:
            [bgzipped_file]
    index:
        run: index.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
    variants_to_table:
        run: variants_to_table.cwl
        in:
            reference: reference
            vcf: index/indexed_vcf
            fields: variants_to_table_fields
            genotype_fields: variants_to_table_genotype_fields
        out:
            [variants_tsv]
    add_vep_fields_to_table:
        run: add_vep_fields_to_table.cwl
        in:
            vcf: index/indexed_vcf
            vep_fields: vep_to_table_fields
            tsv: variants_to_table/variants_tsv
        out: [annotated_variants_tsv]
