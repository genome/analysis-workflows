#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "CLE gold vcf evaluation workflow"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    roi_bed:
        type: File
    query_vcf:
        type: File
        secondaryFiles: [.tbi]
    normal_cram:
        type: File
        secondaryFiles: [^.crai]
    tumor_cram:
        type: File
        secondaryFiles: [^.crai]
    normal_sample_name:
        type: string?
        default: 'NORMAL'
    tumor_sample_name:
        type: string?
        default: 'TUMOR'
    gold_vcf:
        type: File
        secondaryFiles: [.tbi]
    true_negative_bed:
        type: File
    output_dir:
        type: string
outputs:
    output_directory:
        type: Directory
        outputSource: gather_to_sub_directory/gathered_directory
steps:
    normal_cram_to_bam_and_index:
        run: cram_to_bam_and_index.cwl
        in:
            cram: normal_cram
            reference: reference
        out:
            [bam]
    tumor_cram_to_bam_and_index:
        run: cram_to_bam_and_index.cwl
        in:
            cram: tumor_cram
            reference: reference
        out:
            [bam]
    gold_vcf_roi:
        run: ../tools/bedtools_intersect.cwl
        in:
            file_a: gold_vcf
            file_b: roi_bed
            output_file_a:
                default: false
            unique_result:
                default: false
            output_name:
                default: 'gold_roi.vcf'
        out:
            [intersect_result]
    bgzip_and_index_gold_roi:
        run: bgzip_and_index.cwl
        in:
            vcf: gold_vcf_roi/intersect_result
        out:
            [indexed_vcf]
    query_vcf_pass:
        run: ../tools/select_variants.cwl
        in:
            reference: reference
            vcf: query_vcf
            exclude_filtered: 
                default: true
        out:
            [filtered_vcf]
    query_vcf_pass_roi:
        run: ../tools/bedtools_intersect.cwl
        in:
            file_a: query_vcf_pass/filtered_vcf
            file_b: roi_bed
            output_file_a:
                default: false
            unique_result:
                default: false
            output_name:
                default: 'query_pass_roi.vcf'
        out:
            [intersect_result]
    bgzip_and_index_query_pass_roi:
        run: bgzip_and_index.cwl
        in:
            vcf: query_vcf_pass_roi/intersect_result
        out:
            [indexed_vcf]
    combine_vcf:
        run: ../tools/combine_variants_concordance.cwl
        in:
            reference: reference
            base_vcf: bgzip_and_index_gold_roi/indexed_vcf
            query_vcf: bgzip_and_index_query_pass_roi/indexed_vcf
        out:
           [combined_vcf]
    normal_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: combine_vcf/combined_vcf
            sample: normal_sample_name
            reference_fasta: reference
            bam: normal_cram_to_bam_and_index/bam
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    tumor_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: combine_vcf/combined_vcf
            sample: tumor_sample_name
            reference_fasta: reference
            bam: tumor_cram_to_bam_and_index/bam
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    query_vcf_pass_snv:
        run: ../tools/select_variants.cwl
        in:
            reference: reference
            vcf: bgzip_and_index_query_pass_roi/indexed_vcf
            select_type: 
                default: 'SNP'
        out:
            [filtered_vcf]
    query_vcf_pass_indel:
        run: ../tools/select_variants.cwl
        in:
            reference: reference
            vcf: bgzip_and_index_query_pass_roi/indexed_vcf
            select_type: 
                default: 'INDEL'
        out:
            [filtered_vcf]
    true_negative_roi:
        run: ../tools/bedtools_intersect.cwl
        in:
            file_a: true_negative_bed
            file_b: roi_bed
            output_file_a:
                default: false
            unique_result:
                default: false
            output_name:
                default: 'tn_roi.bed'
        out:
            [intersect_result]
    true_negative_intersect_query_snv:
        run: ../tools/bedtools_intersect.cwl
        in:
            file_a: true_negative_roi/intersect_result
            file_b: query_vcf_pass_snv/filtered_vcf
            output_file_a:
                default: false
            unique_result:
                default: false
            output_name:
                default: 'tn_x_query_snv.bed'
        out:
            [intersect_result]
    true_negative_intersect_query_indel:
        run: ../tools/bedtools_intersect.cwl
        in:
            file_a: true_negative_roi/intersect_result
            file_b: query_vcf_pass_indel/filtered_vcf
            output_file_a:
                default: false
            unique_result:
                default: false
            output_name: 
                default: 'tn_x_query_indel.bed'
        out:
            [intersect_result]
    sompy:
        run: ../tools/sompy.cwl
        in:
            reference: reference
            roi_bed:  roi_bed
            truth_vcf: gold_vcf
            query_vcf: query_vcf_pass/filtered_vcf
        out:
            [sompy_out]
    evaluation:
        run: ../tools/eval_cle_gold.cwl
        in:
            sompy_out: sompy/sompy_out
            true_negative_bed: true_negative_roi/intersect_result
            tn_x_query_snv: true_negative_intersect_query_snv/intersect_result
            tn_x_query_indel: true_negative_intersect_query_indel/intersect_result
            combined_vcf: combine_vcf/combined_vcf
            normal_snv_bam_readcount_tsv: normal_bam_readcount/snv_bam_readcount_tsv
            normal_indel_bam_readcount_tsv: normal_bam_readcount/indel_bam_readcount_tsv
            tumor_snv_bam_readcount_tsv: tumor_bam_readcount/snv_bam_readcount_tsv
            tumor_indel_bam_readcount_tsv: tumor_bam_readcount/indel_bam_readcount_tsv
        out:
            [eval_out, vaf_report]
    gather_to_sub_directory:
        run: ../tools/gather_to_sub_directory.cwl
        in:
            outdir: output_dir
            files: [evaluation/eval_out,evaluation/vaf_report]
        out:
            [gathered_directory]
