#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: Workflow
label: "Per-interval mutect"
inputs:
    reference:
        type: File
        secondaryFiles: [".fai", "^.dict"]
    tumor_bam:
        type: File
        secondaryFiles: .bai
    normal_bam:
        type: File
        secondaryFiles: .bai
    interval_list:
        type: File
outputs:
    vcf:
        type: File
        outputSource: index/indexed_vcf
steps:
    mutect:
        run: mutect.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: interval_list
        out:
            [vcf]
    index:
        run: index.cwl
        in:
            vcf: [mutect/vcf]
        out:
            [indexed_vcf]

