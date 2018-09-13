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
outputs:
    unfiltered_vcf:
        type: File
        outputSource: GATK_haplotype_caller/docm_out
    filtered_vcf:
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
            [docm_out]
    docm_filter:
        run: ../tools/somatic_docm_filter.cwl
        in:
            docm_out: GATK_haplotype_caller/docm_out
            normal_cram: normal_cram
            tumor_cram: tumor_cram
        out:
            [docm_filter_out]
    bgzip:
        run: ../tools/bgzip.cwl
        in:
            file: docm_filter/docm_filter_out
        out:
            [bgzipped_file]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
