#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Get snv and indel readcounts for one sample and one bam and add them to a vcf"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    vcf:
        type: File
        secondaryFiles: [.tbi]
    sample_name:
        type: string
    reference_fasta:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    bam:
        type: File
        secondaryFiles: [.bai, .crai]
    min_base_quality:
        type: int?
    min_mapping_quality:
        type: int?
    data_type:
        type:
            - type: enum
              symbols: ["DNA", "RNA"]
outputs:
    readcount_vcf:
        type: File
        outputSource: add_indel_bam_readcount_to_vcf/annotated_bam_readcount_vcf
steps:
    normalize_variants:
        run: ../tools/normalize_variants.cwl
        in:
            reference: reference_fasta
            vcf: vcf
        out:
            [normalized_vcf]
    decompose_variants:
        run: ../tools/vt_decompose.cwl
        in:
            vcf: normalize_variants/normalized_vcf
        out:
            [decomposed_vcf]
    bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: decompose_variants/decomposed_vcf
            sample: sample_name
            reference_fasta: reference_fasta
            bam: bam
            min_base_quality: min_base_quality
            min_mapping_quality: min_mapping_quality
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    add_snv_bam_readcount_to_vcf:
        run: ../tools/vcf_readcount_annotator.cwl
        in:
            vcf: vcf
            bam_readcount_tsv: bam_readcount/snv_bam_readcount_tsv
            data_type: data_type
            sample_name: sample_name
            variant_type:
                default: 'snv'
        out:
            [annotated_bam_readcount_vcf]
    add_indel_bam_readcount_to_vcf:
        run: ../tools/vcf_readcount_annotator.cwl
        in:
            vcf: add_snv_bam_readcount_to_vcf/annotated_bam_readcount_vcf
            bam_readcount_tsv: indel_bam_readcount_tsv
            data_type: data_type
            sample_name: sample_name
            variant_type:
                default: 'indel'
        out:
            [annotated_bam_readcount_vcf]
