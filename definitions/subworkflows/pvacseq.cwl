#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Workflow to run pVACseq from detect_variants and rnaseq pipeline outputs"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    detect_variants_vcf:
        type: File
    sample_name:
        type: string?
        default: 'TUMOR'
    normal_sample_name:
        type: string?
        default: 'NORMAL'
    rnaseq_bam:
        type: File
        secondaryFiles: ['.bai']
    reference_fasta:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
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
    epitope_lengths_class_i:
        type: int[]?
    epitope_lengths_class_ii:
        type: int[]?
    binding_threshold:
        type: int?
    percentile_threshold:
        type: int?
    allele_specific_binding_thresholds:
        type: boolean?
    minimum_fold_change:
        type: float?
    peptide_sequence_length:
        type: int?
    top_score_metric:
        type:
            - "null"
            - type: enum
              symbols: ["lowest", "median"]
    additional_report_columns:
        type:
            - "null"
            - type: enum
              symbols: ["sample_name"]
    fasta_size:
        type: int?
    downstream_sequence_length:
        type: string?
    exclude_nas:
        type: boolean?
    phased_proximal_variants_vcf:
        type: File?
        secondaryFiles: ['.tbi']
    maximum_transcript_support_level:
        type:
            - "null"
            - type: enum
              symbols: ["1", "2", "3", "4", "5"]
    normal_cov:
        type: int?
    tdna_cov:
        type: int?
    trna_cov:
        type: int?
    normal_vaf:
        type: float?
    tdna_vaf:
        type: float?
    trna_vaf:
        type: float?
    expn_val:
        type: float?
    net_chop_method:
        type:
            - "null"
            - type: enum
              symbols: ["cterm", "20s"]
    net_chop_threshold:
        type: float?
    netmhc_stab:
        type: boolean?
    run_reference_proteome_similarity:
        type: boolean?
    n_threads:
        type: int?
    variants_to_table_fields:
        type: string[]?
        default: [CHROM,POS,ID,REF,ALT]
    variants_to_table_genotype_fields:
        type: string[]?
        default: [GT,AD,AF,DP,RAD,RAF,RDP,GX,TX]
    vep_to_table_fields:
        type: string[]?
        default: [HGVSc,HGVSp]
outputs:
    annotated_vcf:
        type: File
        outputSource: add_transcript_expression_data_to_vcf/annotated_expression_vcf
    annotated_tsv:
        type: File
        outputSource: add_vep_fields_to_table/annotated_variants_tsv
    pvacseq_predictions:
        type: Directory
        outputSource: pvacseq/pvacseq_predictions
steps:
    tumor_rna_bam_readcount:
        run: bam_readcount.cwl
        in:
            vcf: detect_variants_vcf
            sample: sample_name
            reference_fasta: reference_fasta
            bam: rnaseq_bam
            min_base_quality: readcount_minimum_base_quality
            min_mapping_quality: readcount_minimum_mapping_quality
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv, normalized_vcf]
    add_tumor_rna_bam_readcount_to_vcf:
        run: vcf_readcount_annotator.cwl
        in:
            vcf: tumor_rna_bam_readcount/normalized_vcf
            snv_bam_readcount_tsv: tumor_rna_bam_readcount/snv_bam_readcount_tsv
            indel_bam_readcount_tsv: tumor_rna_bam_readcount/indel_bam_readcount_tsv
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
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: add_transcript_expression_data_to_vcf/annotated_expression_vcf
        out:
            [indexed_vcf]
    pvacseq:
        run: ../tools/pvacseq.cwl
        in:
            input_vcf: index/indexed_vcf
            sample_name: sample_name
            alleles: alleles
            prediction_algorithms: prediction_algorithms
            epitope_lengths_class_i: epitope_lengths_class_i
            epitope_lengths_class_ii: epitope_lengths_class_ii
            binding_threshold: binding_threshold
            percentile_threshold: percentile_threshold
            normal_sample_name: normal_sample_name
            minimum_fold_change: minimum_fold_change
            top_score_metric: top_score_metric
            additional_report_columns: additional_report_columns
            fasta_size: fasta_size
            downstream_sequence_length: downstream_sequence_length
            exclude_nas: exclude_nas
            phased_proximal_variants_vcf: phased_proximal_variants_vcf
            maximum_transcript_support_level: maximum_transcript_support_level
            normal_cov: normal_cov
            tdna_cov: tdna_cov
            trna_cov: trna_cov
            normal_vaf: normal_vaf
            tdna_vaf: tdna_vaf
            trna_vaf: trna_vaf
            expn_val: expn_val
            net_chop_method: net_chop_method
            net_chop_threshold: net_chop_threshold
            netmhc_stab: netmhc_stab
            run_reference_proteome_similarity: run_reference_proteome_similarity
            n_threads: n_threads
        out:
            [pvacseq_predictions]
    variants_to_table:
        run: ../tools/variants_to_table.cwl
        in:
            reference: reference_fasta
            vcf: index/indexed_vcf
            fields: variants_to_table_fields
            genotype_fields: variants_to_table_genotype_fields
        out:
            [variants_tsv]
    add_vep_fields_to_table:
        run: ../tools/add_vep_fields_to_table.cwl
        in:
            vcf: index/indexed_vcf
            vep_fields: vep_to_table_fields
            tsv: variants_to_table/variants_tsv
            prefix:
                default: 'pvacseq'
        out: [annotated_variants_tsv]
