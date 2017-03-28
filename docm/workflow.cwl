#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Detect Docm variants"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type: File
        secondaryFiles: [.fai]
    tumor_cram:
        type: File
        secondaryFiles: [^.crai]
    normal_cram:
        type: File
        secondaryFiles: [^.crai]
    docm_vcf:
        type: File
        secondaryFiles: [.tbi]
outputs:
    merged_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [.tbi]
steps:
    GATK_haplotype_caller:
        run: GATK_haplotype_caller.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            docm_vcf: docm_vcf
            interval_file: docm_vcf
        out:
            [docm_out]
    docm_filter:
        run: docm_filter.cwl
        in:
            docm_out: GATK_haplotype_caller/docm_out
        out:
            [docm_filter_out]
    bgzip:
        run: ../detect_variants/bgzip.cwl
        in:
            file: docm_filter/docm_filter_out
        out:
            [bgzipped_file]
    index:
        run: ../detect_variants/index.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
