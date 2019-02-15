#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Detect Docm variants"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type: string
    tumor_cram:
        type: File
        secondaryFiles: [^.crai]
    normal_cram:
        type: File
        secondaryFiles: [^.crai]
    docm_vcf:
        type: File
        secondaryFiles: [.tbi]
    interval_list:
        type: File
    filter_docm_variants:
        type: Integer
outputs:
    docm_variants_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [.tbi]
steps:
    GATK_haplotype_caller:
        run: ../tools/docm_gatk_haplotype_caller.cwl
        in:
            reference: reference
            cram: tumor_cram
            normal_cram: normal_cram
            docm_vcf: docm_vcf
            interval_list: interval_list
        out:
            [docm_raw_variants]
    decompose:
        run: ../tools/vt_decompose.cwl
        in:
            vcf: GATK_haplotype_caller/docm_raw_variants
        out:
            [decomposed_vcf]
    docm_filter:
        run: ../tools/somatic_docm_filter.cwl
        in:
            docm_raw_variants: decompose/decomposed_vcf
            normal_cram: normal_cram
            tumor_cram: tumor_cram
            filter_docm_variants: filter_docm_variants
        out:
            [docm_filtered_variants]
    bgzip:
        run: ../tools/bgzip.cwl
        in:
            file: docm_filter/docm_filtered_variants
        out:
            [bgzipped_file]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
