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
    reference: 
        type: string
outputs: 
    filtered_vcf:
        type: File
        outputSource: filter_vcf_cle/cle_filtered_vcf
        secondaryFiles: [.tbi]
steps:
    filter_vcf_gnomADe_allele_freq:
        run: filter_vcf_gnomADe_allele_freq.cwl
        in:
            vcf: vcf
            maximum_population_allele_frequency: filter_gnomADe_maximum_population_allele_frequency
        out:
            [filtered_vcf]
    filter_vcf_mapq0:
        run: filter_vcf_mapq0.cwl
        in: 
            vcf: filter_vcf_gnomADe_allele_freq/filtered_vcf
            tumor_bam: tumor_bam
            threshold: filter_mapq0_threshold
            reference: reference
        out:
            [mapq0_filtered_vcf]
    filter_vcf_cle:
        run: ../definitions/tools/filter_vcf_cle.cwl
        in:
            vcf: filter_vcf_mapq0/mapq0_filtered_vcf
            filter: do_cle_vcf_filter
        out:
            [cle_filtered_vcf]
