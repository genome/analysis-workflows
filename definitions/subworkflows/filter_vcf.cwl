#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Apply filters to VCF file"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: MultipleInputFeatureRequirement
inputs:
    vcf:
        type: File
        secondaryFiles: [.tbi]
    filter_mapq0_threshold:
        type: float
    filter_gnomADe_maximum_population_allele_frequency:
        type: float
    gnomad_field_name:
        type: string
    tumor_bam:
        type: File
        secondaryFiles: [.bai]
    do_cle_vcf_filter:
        type: boolean
    filter_somatic_llr_threshold:
        type: float
    filter_somatic_llr_tumor_purity:
        type: float
    filter_somatic_llr_normal_contamination_rate:
        type: float
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    filter_minimum_depth:
        type: int
    tumor_sample_name:
        type: string
    normal_sample_name:
        type: string
    known_variants:
        type: File?
        secondaryFiles: [.tbi]
        doc: "Previously discovered variants to be flagged in this workflow's output vcf"
outputs:
    filtered_vcf:
        type: File
        outputSource: set_final_vcf_name/replacement
steps:
    filter_known_variants:
        run: ../tools/filter_known_variants.cwl
        in:
            known_variants: known_variants
            vcf: vcf
        out:
            [known_filtered]
    filter_vcf_gnomADe_allele_freq:
        run: ../tools/filter_vcf_custom_allele_freq.cwl
        in:
            vcf: filter_known_variants/known_filtered
            maximum_population_allele_frequency: filter_gnomADe_maximum_population_allele_frequency
            field_name: gnomad_field_name
        out:
            [filtered_vcf]
    filter_vcf_mapq0:
        run: ../tools/filter_vcf_mapq0.cwl
        in:
            vcf: filter_vcf_gnomADe_allele_freq/filtered_vcf
            tumor_bam: tumor_bam
            threshold: filter_mapq0_threshold
            reference: reference
        out:
            [mapq0_filtered_vcf]
    filter_vcf_cle:
        run: ../tools/filter_vcf_cle.cwl
        in:
            vcf: filter_vcf_mapq0/mapq0_filtered_vcf
            filter: do_cle_vcf_filter
        out:
            [cle_filtered_vcf]
    filter_vcf_depth:
        run: ../tools/filter_vcf_depth.cwl
        in:
            vcf: filter_vcf_cle/cle_filtered_vcf
            minimum_depth: filter_minimum_depth
            sample_names:
                source: [normal_sample_name, tumor_sample_name]
                linkMerge: merge_flattened
        out:
            [depth_filtered_vcf]
    filter_vcf_somatic_llr:
        run: ../tools/filter_vcf_somatic_llr.cwl
        in:
            vcf: filter_vcf_depth/depth_filtered_vcf
            threshold: filter_somatic_llr_threshold
            tumor_purity: filter_somatic_llr_tumor_purity
            normal_contamination_rate: filter_somatic_llr_normal_contamination_rate
            tumor_sample_name: tumor_sample_name
            normal_sample_name: normal_sample_name
        out:
            [somatic_llr_filtered_vcf]
    set_final_vcf_name:
        run: ../tools/staged_rename.cwl
        in:
            original: filter_vcf_somatic_llr/somatic_llr_filtered_vcf
            name:
                valueFrom: 'annotated_filtered.vcf'
        out:
            [replacement]
