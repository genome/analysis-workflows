#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "fp_filter workflow"
inputs:
    cram:
        type: File
        secondaryFiles: [^.crai]
    reference:
        type: File
        secondaryFiles: [.fai]
    vcf:
        type: File
        secondaryFiles: [.tbi]
outputs:
    filtered_vcf:
        type: File
        outputSource: fp_filter/filtered_vcf
steps:
    cram_to_bam:
        run: cram_to_bam.cwl
        in:
            cram: cram
            reference: reference
        out:
            [bam]
    fp_filter:
        run: fp_filter.cwl
        in:
            reference: reference
            bam: cram_to_bam/bam
            vcf: vcf
        out:
            [filtered_vcf]


