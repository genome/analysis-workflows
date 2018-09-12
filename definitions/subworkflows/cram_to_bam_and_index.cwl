#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "cram_to_bam workflow"
inputs:
    cram:
        type: File
        secondaryFiles: [^.crai]
    reference:
        type: string
outputs:
    bam:
        type: File
        outputSource: index_bam/indexed_bam
        secondaryFiles: [.bai]
steps:
    cram_to_bam:
        run: ../tools/cram_to_bam.cwl
        in:
            cram: cram
            reference: reference
        out:
            [bam]
    index_bam:
        run: ../tools/index_bam.cwl
        in:
            bam: cram_to_bam/bam
        out:
            [indexed_bam]

