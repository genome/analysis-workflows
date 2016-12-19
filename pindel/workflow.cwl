#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "pindel workflow"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type: File
        secondaryFiles: [".fai"]
    tumor_bam:
        type: File
        secondaryFiles: .bai
    normal_bam:
        type: File
        secondaryFiles: .bai
    interval_list:
        type: File
outputs:
    merged_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: .tbi
steps:
    pindel:
        run: pindel.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: interval_list
        out:
            [deletions, insertions]
    cat:
        run: cat.cwl
        in:
            deletion_out: pindel/deletions
            insertion_out: pindel/insertions
        out:
            [pindel_out]
    grep:
        run: grep.cwl
        in: 
            pindel_output: cat/pindel_out
        out:
            [pindel_output_summary]
    somaticfilter:
        run: somaticfilter.cwl
        in:
            reference: reference
            pindel_output_summary: grep/pindel_output_summary
        out: 
            [vcf]
    bgzip:
        run: ../bgzip.cwl
        in: 
            file: somaticfilter/vcf
        out:
            [bgzipped_file]
    index:
        run: ../index.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
