#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: 'Unaligned to aligned BAM'
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement
inputs:
    bam:
        type: File
    reference:
        type: string
    readgroup:
        type: string
outputs:
    tagged_bam:
        type: File
        outputSource: align_and_tag/aligned_bam
steps:
    align_and_tag:
        run: ../tools/align_and_tag.cwl
        in:
            bam: bam
            readgroup: readgroup
            reference: reference
        out:
            [aligned_bam]
