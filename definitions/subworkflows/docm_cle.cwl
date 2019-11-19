#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Detect Docm variants"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    tumor_bam:
        type: File
        secondaryFiles: [^.bai]
    normal_bam:
        type: File
        secondaryFiles: [^.bai]
    docm_vcf:
        type: File
        secondaryFiles: [.tbi]
    interval_list:
        type: File
    filter_docm_variants:
        type: boolean
outputs:
    docm_variants_vcf:
        type: File
        outputSource: index2/indexed_vcf
        secondaryFiles: [.tbi]
steps:
    GATK_haplotype_caller:
        run: ../tools/docm_gatk_haplotype_caller.cwl
        in:
            reference: reference
            bam: tumor_bam
            normal_bam: normal_bam
            docm_vcf: docm_vcf
            interval_list: interval_list
        out:
            [docm_raw_variants]
    bgzip:
        run: ../tools/bgzip.cwl
        in:
            file: GATK_haplotype_caller/docm_raw_variants
        out:
            [bgzipped_file]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
    decompose:
        run: ../tools/vt_decompose.cwl
        in:
            vcf: index/indexed_vcf
        out:
            [decomposed_vcf]
    docm_filter:
        run: ../tools/filter_vcf_docm.cwl
        in:
            docm_raw_variants: decompose/decomposed_vcf
            normal_bam: normal_bam
            tumor_bam: tumor_bam
            filter_docm_variants: filter_docm_variants
        out:
            [docm_filtered_variants]
    bgzip2:
        run: ../tools/bgzip.cwl
        in:
            file: docm_filter/docm_filtered_variants
        out:
            [bgzipped_file]
    index2:
        run: ../tools/index_vcf.cwl
        in:
            vcf: bgzip2/bgzipped_file
        out:
            [indexed_vcf]
