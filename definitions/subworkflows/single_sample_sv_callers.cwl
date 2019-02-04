#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Subworkflow to allow calling different SV callers which require bam files as inputs"

requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
inputs:
    cram:
        type: File
    reference:
        type: string

    cnvkit_diagram:
        type: boolean?
    cnvkit_drop_low_coverage:
        type: boolean?
    cnvkit_method:
        type: string?
    cnvkit_reference_cnn:
        type: File
    cnvkit_scatter_plot:
        type: boolean?
    cnvkit_male_reference:
        type: boolean?
    cnvkit_vcf_name:
        type: string?

    manta_call_regions:
        type: File?
    manta_non_wgs:
        type: boolean?
    manta_output_contigs:
        type: boolean?

    smoove_exclude_regions:
        type: File?

outputs:
    cn_diagram:
        type: File?
        outputSource: run_cnvkit/cn_diagram
    cn_scatter_plot:
        type: File?
        outputSource: run_cnvkit/cn_scatter_plot
    tumor_antitarget_coverage:
        type: File
        outputSource: run_cnvkit/tumor_antitarget_coverage
    tumor_target_coverage:
        type: File
        outputSource: run_cnvkit/tumor_target_coverage
    tumor_bin_level_ratios:
        type: File
        outputSource: run_cnvkit/tumor_bin_level_ratios
    tumor_segmented_ratios:
        type: File
        outputSource: run_cnvkit/tumor_segmented_ratios
    cnvkit_vcf:
        type: File
        outputSource: run_cnvkit/cnvkit_vcf
    manta_diploid_variants:
        type: File?
        outputSource: run_manta/diploid_variants
    manta_somatic_variants:
        type: File?
        outputSource: run_manta/somatic_variants
    manta_all_candidates:
        type: File
        outputSource: run_manta/all_candidates
    manta_small_candidates:
        type: File
        outputSource: run_manta/small_candidates
    manta_tumor_only_variants:
        type: File?
        outputSource: run_manta/tumor_only_variants
    smoove_output_variants:
        type: File
        outputSource: run_smoove/output_vcf

steps:
    cram_to_bam:
        run: ../tools/cram_to_bam.cwl
        in:
            cram: cram
            reference: reference
        out:
            [bam]
    index_bam:
        run: ../tools/index_bam.cwl
        in:
            bam: cram_to_bam/bam
        out:
            [indexed_bam]
    run_cnvkit:
        run: cnvkit_single_sample.cwl
        in:
            diagram: cnvkit_diagram
            drop_low_coverage: cnvkit_drop_low_coverage
            method: cnvkit_method
            reference_cnn: cnvkit_reference_cnn
            tumor_bam: index_bam/indexed_bam
            scatter_plot: cnvkit_scatter_plot
            male_reference: cnvkit_male_reference
            cnvkit_vcf_name: cnvkit_vcf_name
        out:
            [cn_diagram, cn_scatter_plot, tumor_antitarget_coverage, tumor_target_coverage, tumor_bin_level_ratios, tumor_segmented_ratios, cnvkit_vcf]
    run_manta:
        run: ../tools/manta_somatic.cwl
        in:
            call_regions: manta_call_regions
            non_wgs: manta_non_wgs
            output_contigs: manta_output_contigs
            reference: reference
            tumor_bam: index_bam/indexed_bam
        out: 
            [diploid_variants, somatic_variants, all_candidates, small_candidates, tumor_only_variants]
    run_smoove:
        run: ../tools/smoove.cwl
        in:
            bams:
                source: [index_bam/indexed_bam]
                linkMerge: merge_flattened
            exclude_regions: smoove_exclude_regions
            reference: reference
        out:
            [output_vcf]

