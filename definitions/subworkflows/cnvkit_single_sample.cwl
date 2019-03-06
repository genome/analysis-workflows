#! /usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Subworkflow that runs cnvkit in single sample mode and returns a vcf file"

inputs:
    tumor_bam:
        type: File
    method:
        type: string?
        default: "hybrid"
    diagram:
        type: boolean?
    scatter_plot:
        type: boolean?
    drop_low_coverage:
        type: boolean?
    male_reference:
        type: boolean?
    reference_cnn:
        type: File
    cnvkit_vcf_name:
        type: string?
        default: $(inputs.tumor_bam.nameroot).cnvkit.vcf
outputs:
    cn_diagram:
        type: File?
        outputSource: cnvkit_main/cn_diagram
    cn_scatter_plot:
        type: File?
        outputSource: cnvkit_main/cn_scatter_plot
    tumor_antitarget_coverage:
        type: File
        outputSource: cnvkit_main/tumor_antitarget_coverage
    tumor_target_coverage:
        type: File
        outputSource: cnvkit_main/tumor_target_coverage
    tumor_bin_level_ratios:
        type: File
        outputSource: cnvkit_main/tumor_bin_level_ratios
    tumor_segmented_ratios:
        type: File
        outputSource: cnvkit_main/tumor_segmented_ratios
    cnvkit_vcf:
        type: File
        outputSource: cns_to_vcf/cnvkit_vcf

steps:
    cnvkit_main:
        run: ../tools/cnvkit_batch.cwl
        in: 
            tumor_bam: tumor_bam
            method: method
            diagram: diagram
            scatter_plot: scatter_plot
            drop_low_coverage: drop_low_coverage
            male_reference: male_reference
            reference_cnn: reference_cnn
        out:
            [cn_diagram, cn_scatter_plot, tumor_antitarget_coverage, tumor_target_coverage, tumor_bin_level_ratios, tumor_segmented_ratios]            
            
    cns_to_vcf:
        run: ../tools/cnvkit_vcf_export.cwl
        in:
            cns_file: cnvkit_main/tumor_segmented_ratios
            male_reference: male_reference
            cnr_file: cnvkit_main/tumor_bin_level_ratios
            output_name: cnvkit_vcf_name
        out:
            [cnvkit_vcf]
