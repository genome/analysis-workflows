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
outputs:
    annotated_vcf:
        type: File
        outputSource: add_tumor_rna_bam_readcount_to_vcf/annotated_bam_readcount_vcf
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
