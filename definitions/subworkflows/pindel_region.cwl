#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Run pindel on provided region"
requirements:
    - class: ScatterFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    reference:
        type: string
    tumor_bam:
        type: File
        secondaryFiles: ["^.bai"]
    normal_bam:
        type: File
        secondaryFiles: ["^.bai"]
    region_file:
        type: File
    insert_size:
        type: int
        default: 400
    ref_name:
        type: string?
        default: "GRCh38DH"
    output_name:
        type: string?
        default: "pindel.vcf"
    tumor_sample_name:
        type: string
    normal_sample_name:
        type: string
outputs:
    pindel_region_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [".tbi"]
steps:
    pindel_region:
        run: ../tools/pindel.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            insert_size: insert_size
            region_file: region_file
            tumor_sample_name: tumor_sample_name
            normal_sample_name: normal_sample_name
        out:
            [deletions, insertions, tandems, long_insertions, inversions]
    cat:
        run: ../tools/cat_out.cwl
        in:
            pindel_outs: [pindel_region/deletions, pindel_region/insertions, pindel_region/tandems, pindel_region/long_insertions, pindel_region/inversions]
        out:
            [pindel_out]
    pindel2vcf:
        run: ../tools/pindel2vcf.cwl
        in:
            pindel_out: cat/pindel_out
            reference: reference
            ref_name:  ref_name
            output_name: output_name
        out:
            [pindel_vcf]
    fix_vcf_header:
        run: ../tools/fix_vcf_header.cwl
        in: 
            vcf: pindel2vcf/pindel_vcf
        out:
            [fixed_vcf]
    remove_end_tags:
        run: ../tools/remove_end_tags.cwl
        in:
            vcf: fix_vcf_header/fixed_vcf
        out:
            [processed_vcf]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: remove_end_tags/processed_vcf
        out:
            [indexed_vcf]
