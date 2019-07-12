#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Apply filters to VCF file"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    vcf:
        type: File
    filter_mapq0_threshold: 
        type: float
    filter_gnomADe_maximum_population_allele_frequency:
        type: float
    tumor_bam: 
        type: File
        secondaryFiles: [.bai]
    do_cle_vcf_filter: 
        type: boolean
    filter_somatic_llr_threshold:
        type: float
    reference: 
        type: string
    filter_minimum_depth:
        type: int
    sample_names:
        type: string
outputs: 
    filtered_vcf:
        type: File
        outputSource: filter_vcf_somatic_llr/somatic_llr_filtered_vcf
steps:
    filter_vcf_gnomADe_allele_freq:
        run: ../tools/filter_vcf_gnomADe_allele_freq.cwl
        in:
            vcf: vcf
            maximum_population_allele_frequency: filter_gnomADe_maximum_population_allele_frequency
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
            sample_names: sample_names
        out:
            [depth_filtered_vcf]
    filter_vcf_somatic_llr:
        run: ../tools/filter_vcf_somatic_llr.cwl
        in:
            vcf: filter_vcf_depth/depth_filtered_vcf
            threshold: filter_somatic_llr_threshold
        out:
            [somatic_llr_filtered_vcf]
