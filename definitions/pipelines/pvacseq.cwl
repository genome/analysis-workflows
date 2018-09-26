#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Workflow to run pVACseq from detect_variants and rnaseq pipeline outputs"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    detect_variants_vcf:
        type: File
        secondaryFiles: ['.tbi']
    sample_name:
        type: string?
        default: 'TUMOR'
    rnaseq_bam:
        type: File
        secondaryFiles: ['.bai']
    reference_fasta:
        type: string
    readcount_minimum_base_quality:
        type: int?
    readcount_minimum_mapping_quality:
        type: int?
    gene_expression_file:
        type: File
    transcript_expression_file:
        type: File
    expression_tool:
        type: string?
        default: 'kallisto'
    alleles:
        type: string[]
    prediction_algorithms:
        type: string[]
    epitope_lengths:
        type: int[]?
outputs:
    mhc_i_all_epitopes:
        type: File?
        outputSource: pvacseq/mhc_i_all_epitopes
    mhc_i_filtered_epitopes:
        type: File?
        outputSource: pvacseq/mhc_i_filtered_epitopes
    mhc_i_ranked_epitopes:
        type: File?
        outputSource: pvacseq/mhc_i_ranked_epitopes
    mhc_ii_all_epitopes:
        type: File?
        outputSource: pvacseq/mhc_ii_all_epitopes
    mhc_ii_filtered_epitopes:
        type: File?
        outputSource: pvacseq/mhc_ii_filtered_epitopes
    mhc_ii_ranked_epitopes:
        type: File?
        outputSource: pvacseq/mhc_ii_ranked_epitopes
    combined_all_epitopes:
        type: File?
        outputSource: pvacseq/combined_all_epitopes
    combined_filtered_epitopes:
        type: File?
        outputSource: pvacseq/combined_filtered_epitopes
    combined_ranked_epitopes:
        type: File?
        outputSource: pvacseq/combined_ranked_epitopes
steps:
    tumor_rna_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: detect_variants_vcf
            sample: sample_name
            reference_fasta: reference_fasta
            bam: rnaseq_bam
            min_base_quality: readcount_minimum_base_quality
            min_mapping_quality: readcount_minimum_mapping_quality
        out:
            [bam_readcount_tsv]
    add_tumor_rna_bam_readcount_to_vcf:
        run: ../tools/vcf_readcount_annotator.cwl
        in:
            vcf: detect_variants_vcf
            bam_readcount_tsv: tumor_rna_bam_readcount/bam_readcount_tsv
            data_type:
                default: 'RNA'
            sample_name: sample_name
        out:
            [annotated_bam_readcount_vcf]
    add_gene_expression_data_to_vcf:
        run: ../tools/vcf_expression_annotator.cwl
        in:
            vcf: add_tumor_rna_bam_readcount_to_vcf/annotated_bam_readcount_vcf
            expression_file: gene_expression_file
            expression_tool: expression_tool
            data_type:
                default: 'gene'
            sample_name: sample_name
        out:
            [annotated_expression_vcf]
    add_transcript_expression_data_to_vcf:
        run: ../tools/vcf_expression_annotator.cwl
        in:
            vcf: add_gene_expression_data_to_vcf/annotated_expression_vcf
            expression_file: transcript_expression_file
            expression_tool: expression_tool
            data_type:
                default: 'transcript'
            sample_name: sample_name
        out:
            [annotated_expression_vcf]
    pvacseq:
        run: ../tools/pvacseq.cwl
        in:
            input_file: add_transcript_expression_data_to_vcf/annotated_expression_vcf
            sample_name: sample_name
            alleles: alleles
            prediction_algorithms: prediction_algorithms
            epitope_lengths: epitope_lengths
        out:
            [
                mhc_i_all_epitopes,
                mhc_i_filtered_epitopes,
                mhc_i_ranked_epitopes,
                mhc_ii_all_epitopes,
                mhc_ii_filtered_epitopes,
                mhc_ii_ranked_epitopes,
                combined_all_epitopes,
                combined_filtered_epitopes,
                combined_ranked_epitopes,
            ]
