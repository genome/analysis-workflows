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
        secondaryFiles: [.tbi]
steps:
    add_gt:
        run: ../definitions/tools/add_strelka_gt.cwl
        in:
            vcf: vcf
        out:
            [processed_vcf]
    bgzip:
        run: ../definitions/tools/bgzip.cwl
        in:
            file: add_gt/processed_vcf
        out:
            [bgzipped_file]
    index:
        run: ../definitions/tools/index_vcf.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
