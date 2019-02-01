#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Add snv and indel bam-readcount files to a vcf"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    vcf:
        type: File
        secondaryFiles: [.tbi]
    snv_bam_readcount_tsv:
        type: File
    indel_bam_readcount_tsv:
        type: File
    data_type:
        type:
            - type: enum
              symbols: ["DNA", "RNA"]
    sample_name:
        type: string
outputs:
    annotated_bam_readcount_vcf:
        type: File
        outputSource: add_indel_bam_readcount_to_vcf/annotated_bam_readcount_vcf
steps:
    add_snv_bam_readcount_to_vcf:
        run: ../tools/vcf_readcount_annotator.cwl
        in:
            vcf: vcf
            bam_readcount_tsv: snv_bam_readcount_tsv
            data_type: data_type
            sample_name: sample_name
            variant_type:
                default: 'snv'
        out:
            [annotated_bam_readcount_vcf]
    add_indel_bam_readcount_to_vcf:
        run: ../tools/vcf_readcount_annotator.cwl
        in:
            vcf: add_snv_bam_readcount_to_vcf/annotated_bam_readcount_vcf
            bam_readcount_tsv: indel_bam_readcount_tsv
            data_type: data_type
            sample_name: sample_name
            variant_type:
                default: 'indel'
        out:
            [annotated_bam_readcount_vcf]
