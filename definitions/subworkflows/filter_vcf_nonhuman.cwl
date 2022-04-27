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
    filter_mapq0_threshold: 
        type: float
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
    filter_minimum_depth:
        type: int
    tumor_sample_name:
        type: string
    normal_sample_name:
        type: string
outputs: 
    filtered_vcf:
        type: File
        outputSource: filter_vcf_somatic_llr/somatic_llr_filtered_vcf
steps:
    filter_vcf_mapq0:
        run: ../tools/filter_vcf_mapq0.cwl
        in: 
            vcf: vcf
            tumor_bam: tumor_bam
            threshold: filter_mapq0_threshold
            sample_name: tumor_sample_name
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
