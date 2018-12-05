#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "phase VCF"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    somatic_vcf:
        type: File
    germline_vcf:
        type: File
    reference:
        type: string
    reference_dict:
        type: File
    bam:
        type: File
outputs:
    phased_vcf:
        type: File
        outputSource: bgzip_and_index/indexed_vcf
        secondaryFiles: [.tbi]
steps:
    combine_variants:
        run: ../tools/pvacseq_combine_variants.cwl
        in:
            reference: reference
            germline_vcf: germline_vcf
            somatic_vcf: somatic_vcf
        out:
            [combined_vcf]
    sort:
        run: ../tools/sort_vcf.cwl
        in:
            vcf: combine_variants/combined_vcf
            reference_dict: reference_dict
        out:
            [sorted_vcf]
    phase_vcf:
        run: ../tools/read_backed_phasing.cwl
        in:
            reference: reference
            bam: bam
            vcf: sort/sorted_vcf
        out:
            [phased_vcf]
    bgzip_and_index:
        run: bgzip_and_index.cwl
        in:
            vcf: phase_vcf/phased_vcf
        out:
            [indexed_vcf]
