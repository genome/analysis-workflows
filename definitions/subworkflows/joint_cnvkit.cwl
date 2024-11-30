#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "jointly run cnvkit for sv calls"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement
    - class: ScatterFeatureRequirement
inputs:
    sample_names:
        type: string[]
    bams:
        type: File[]
        secondaryFiles: [^.bai]
    reference_fasta:
        type:
            - string
            - File
        secondaryFiles: [.fai]
    reference_cnn:
        type: File?
        doc: "can be a flat reference or reference based on a panel of normals"
    method:
        type:
          - "null"
          - type: enum
            symbols: ["hybrid", "amplicon", "wgs"]
    segment_filter:
        type:
          - "null"
          - type: enum
            symbols: ["ampdel", "ci", "cn", "sem"]
outputs:
    vcfs:
        type: File[]
        outputSource: index_cnvkit/indexed_vcf
        secondaryFiles: [.tbi]
    cnr:
        type: File[]
        outputSource: cnvkit/tumor_bin_level_ratios
    cns:
        type: File[]
        outputSource: cnvkit/tumor_segmented_ratios
steps:
    cnvkit:
        scatter: [tumor_bam, cnvkit_vcf_name]
        scatterMethod: dotproduct
        run: cnvkit_single_sample.cwl
        in:
            method: method
            reference_cnn: reference_cnn
            tumor_bam: bams
            cnvkit_vcf_name:
                source: [sample_names]
                valueFrom: "$(self).cnvkit.vcf"
            segment_filter: segment_filter
            fasta_reference: reference_fasta
        out:
            [tumor_bin_level_ratios, tumor_segmented_ratios, cnvkit_vcf]
    bgzip_and_index:
        scatter: [vcf]
        run: bgzip_and_index.cwl
        in:
            vcf: cnvkit/cnvkit_vcf
        out:
            [indexed_vcf]
    sample_rename:
        scatter: [input_vcf, new_sample_name]
        scatterMethod: dotproduct
        run: ../tools/replace_vcf_sample_name.cwl
        in:
            input_vcf: bgzip_and_index/indexed_vcf
            new_sample_name: sample_names
            sample_to_replace:
                valueFrom: 'adjusted.tumor'
            output_name:
                valueFrom: '${
                    var sample = inputs.new_sample_name;
                    var name = sample + ".cnvkit.vcf.gz";
                    return name;
                }'
        out:
            [renamed_vcf]
    index_cnvkit:
        scatter: [vcf]
        run: ../tools/index_vcf.cwl
        in:
            vcf: sample_rename/renamed_vcf
        out:
            [indexed_vcf]
