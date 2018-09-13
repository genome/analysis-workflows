#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Detect DoCM variants"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type: string
    cram:
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
        outputSource: gatk_haplotypecaller/docm_out
    filtered_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [.tbi]
steps:
    gatk_haplotypecaller:
        run: ../tools/docm_gatk_haplotype_caller.cwl
        in:
            reference: reference
            cram: cram
            docm_vcf: docm_vcf
            interval_list: interval_list
        out:
            [docm_out]
    docm_filter:
        run: ../tools/single_sample_docm_filter.cwl
        in:
            docm_out: gatk_haplotypecaller/docm_out
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
