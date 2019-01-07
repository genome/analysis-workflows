#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "bam_readcount workflow"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    vcf:
        type: File
        secondaryFiles: [.tbi]
    sample:
        type: string
    reference_fasta:
        type: string
    bam:
        type: File
        secondaryFiles: [.bai]
    min_base_quality:
        type: int?
    min_mapping_quality:
        type: int?
outputs:
    normalized_vcf:
        type: File
        outputSource: decompose_variants/decomposed_vcf
    snv_bam_readcount_tsv:
        type: File
        outputSource: bam_readcount/snv_bam_readcount_tsv
    indel_bam_readcount_tsv:
        type: File
        outputSource: bam_readcount/indel_bam_readcount_tsv
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
            sample: sample
            reference_fasta: reference_fasta
            bam: bam
            min_base_quality: min_base_quality
            min_mapping_quality: min_mapping_quality
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
