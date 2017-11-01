#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "downsample unaligned BAM and align"
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement

inputs:
    bams:
        type: File[]
    readgroups:
        type: string[]
    reference:
        type: string
    dbsnp:
        type: File
        secondaryFiles: [.tbi]
    mills:
        type: File
        secondaryFiles: [.tbi]
    known_indels:
        type: File
        secondaryFiles: [.tbi]
    downsample_probability:
        type: float
outputs:
    cram:
        type: File
        outputSource: align_workflow/cram
        secondaryFiles: [.crai, ^.crai]
steps:
    downsample:
        scatter: [bam]
        scatterMethod: dotproduct
        run: downsample.cwl
        in:
            bam: bams
            probability: downsample_probability
            reference: reference
        out:
            [downsampled_bam]
    align_workflow:
        run: workflow.cwl
        in:
            bams: downsample/downsampled_bam
            readgroups: readgroups
            reference: reference
            dbsnp: dbsnp
            mills: mills
            known_indel: known_indels
        out:
            [cram]
