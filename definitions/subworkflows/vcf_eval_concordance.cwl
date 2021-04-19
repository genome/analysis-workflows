#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Vcf concordance evaluation workflow"
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
    base_normal_cram:
        type: File
        secondaryFiles: [^.crai]
    base_tumor_cram:
        type: File
        secondaryFiles: [^.crai]
    base_vcf:
        type: File
        secondaryFiles: [.tbi]
    query_normal_cram:
        type: File
        secondaryFiles: [^.crai]
    query_tumor_cram:
        type: File
        secondaryFiles: [^.crai]
    query_vcf:
        type: File
        secondaryFiles: [.tbi]
    normal_sample_name:
        type: string?
        default: 'NORMAL'
    tumor_sample_name:
        type: string?
        default: 'TUMOR'
    output_dir:
        type: string
outputs:
    output_directory:
        type: Directory
        outputSource: gather_to_sub_directory/gathered_directory
steps:
    base_normal_cram_to_bam_and_index:
        run: cram_to_bam_and_index.cwl
        in:
            cram: base_normal_cram
            reference: reference
        out:
            [bam]
    base_tumor_cram_to_bam_and_index:
        run: cram_to_bam_and_index.cwl
        in:
            cram: base_tumor_cram
            reference: reference
        out:
            [bam]
    query_normal_cram_to_bam_and_index:
        run: cram_to_bam_and_index.cwl
        in:
            cram: query_normal_cram
            reference: reference
        out:
            [bam]
    query_tumor_cram_to_bam_and_index:
        run: cram_to_bam_and_index.cwl
        in:
            cram: query_tumor_cram
            reference: reference
        out:
            [bam]
    base_vcf_pass:
        run: ../tools/select_variants.cwl
        in:
            reference: reference
            vcf: base_vcf
            exclude_filtered: 
                default: true
        out:
            [filtered_vcf]
    base_vcf_pass_roi:
        run: ../tools/bedtools_intersect.cwl
        in:
            file_a: base_vcf_pass/filtered_vcf
            file_b: roi_bed
            output_file_a:
                default: false
            unique_result:
                default: false
            output_name:
                default: 'base_pass_roi.vcf'
        out:
            [intersect_result]
    bgzip_and_index_base_pass_roi:
        run: bgzip_and_index.cwl
        in:
            vcf: base_vcf_pass_roi/intersect_result
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
            base_vcf: bgzip_and_index_base_pass_roi/indexed_vcf
            query_vcf: bgzip_and_index_query_pass_roi/indexed_vcf
        out:
            [combined_vcf]
    base_normal_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: combine_vcf/combined_vcf
            sample: normal_sample_name
            reference_fasta: reference
            bam: base_normal_cram_to_bam_and_index/bam
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    base_tumor_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: combine_vcf/combined_vcf
            sample: tumor_sample_name
            reference_fasta: reference
            bam: base_tumor_cram_to_bam_and_index/bam
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    query_normal_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: combine_vcf/combined_vcf
            sample: normal_sample_name
            reference_fasta: reference
            bam: query_normal_cram_to_bam_and_index/bam
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    query_tumor_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: combine_vcf/combined_vcf
            sample: tumor_sample_name
            reference_fasta: reference
            bam: query_tumor_cram_to_bam_and_index/bam
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    sompy:
         run: ../tools/sompy.cwl
         in:
            reference: reference
            roi_bed:  roi_bed
            truth_vcf: base_vcf_pass/filtered_vcf
            query_vcf: query_vcf_pass/filtered_vcf
         out:
            [sompy_out]
    vaf_report:
        run: ../tools/eval_vaf_report.cwl
        in:
            combined_vcf: combine_vcf/combined_vcf
            base_normal_snv_bam_readcount_tsv: base_normal_bam_readcount/snv_bam_readcount_tsv
            base_normal_indel_bam_readcount_tsv: base_normal_bam_readcount/indel_bam_readcount_tsv
            base_tumor_snv_bam_readcount_tsv: base_tumor_bam_readcount/snv_bam_readcount_tsv
            base_tumor_indel_bam_readcount_tsv: base_tumor_bam_readcount/indel_bam_readcount_tsv
            query_normal_snv_bam_readcount_tsv: query_normal_bam_readcount/snv_bam_readcount_tsv
            query_normal_indel_bam_readcount_tsv: query_normal_bam_readcount/indel_bam_readcount_tsv
            query_tumor_snv_bam_readcount_tsv: query_tumor_bam_readcount/snv_bam_readcount_tsv
            query_tumor_indel_bam_readcount_tsv: query_tumor_bam_readcount/indel_bam_readcount_tsv
        out:
            [out_file]
    plot_output:
        run: ../tools/somatic_concordance_graph.cwl
        in:
            sompy_file: sompy/sompy_out
            vaf_file: vaf_report/out_file
        out:
            [out_pdf]
    gather_to_sub_directory:
        run: ../tools/gather_to_sub_directory.cwl
        in:
            outdir: output_dir
            files: [sompy/sompy_out, vaf_report/out_file, plot_output/out_pdf]
        out:
            [gathered_directory]
