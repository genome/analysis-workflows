#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Apply filters to VCF file"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    annotated_vcf:
        type: File
    filter_gnomAD_maximum_population_allele_frequency:
        type: float
    gnomad_field_name:
        type: string
    limit_variant_intervals:
        type: File
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
outputs: 
    filtered_vcf:
        type: File
        secondaryFiles: [.tbi]
        outputSource: index_filtered_vcf/indexed_vcf
    final_vcf:
        type: File
        secondaryFiles: [.tbi]
        outputSource: limit_variants/filtered_vcf
steps:
    coding_variant_filter:
        run: ../tools/filter_vcf_coding_variant.cwl
        in:
            vcf: annotated_vcf
        out:
            [filtered_vcf]
    gnomad_frequency_filter:
        run: ../tools/filter_vcf_custom_allele_freq.cwl
        in:
            vcf: coding_variant_filter/filtered_vcf
            maximum_population_allele_frequency: filter_gnomAD_maximum_population_allele_frequency
            field_name: gnomad_field_name
        out:
            [filtered_vcf]
    set_filtered_vcf_name:
        run: ../tools/staged_rename.cwl
        in:
            original: gnomad_frequency_filter/filtered_vcf
            name:
                valueFrom: 'annotated.filtered.vcf'
        out:
            [replacement]
    bgzip_filtered_vcf:
        run: ../tools/bgzip.cwl
        in:
            file: set_filtered_vcf_name/replacement
        out:
            [bgzipped_file]
    index_filtered_vcf:
        run: ../tools/index_vcf.cwl
        in:
            vcf: bgzip_filtered_vcf/bgzipped_file
        out:
            [indexed_vcf]
    limit_variants:
        run: ../tools/select_variants.cwl
        in:
            reference: reference
            vcf: index_filtered_vcf/indexed_vcf
            interval_list: limit_variant_intervals
            exclude_filtered:
                default: true
            output_vcf_basename:
                default: 'annotated.filtered.final'
        out:
            [filtered_vcf]
