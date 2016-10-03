#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "process VCF workflow"
inputs:
    vcf:
        type: File
outputs:
    processed_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: .tbi
steps:
    add_gt:
        run: add_gt.cwl
        in:
            vcf: vcf
        out:
            [processed_vcf]
    bgzip:
        run: ../detect_variants/bgzip.cwl
        in:
            file: add_gt/processed_vcf
        out:
            [bgzipped_file]
    index:
        run: ../detect_variants/index.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
