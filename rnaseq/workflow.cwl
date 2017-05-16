#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "RnaSeq alignment workflow"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
inputs:
    reference_index:
        type: File #this requires an extra file with the basename
        secondaryFiles: [".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"]
    instrument_data_bams:
        type: File[]
    read_group_id:
        type: string[]
    read_group_fields:
        type:
            type: array
            items:
                type: array
                items: string
outputs:
    aligned_bam:
        type: File
        outputSource: merge/merged_bam
steps:
    align:
        run: align.cwl
        scatter: [instrument_data_bam, read_group_id, read_group_fields]
        scatterMethod: dotproduct
        in:
            instrument_data_bam: instrument_data_bams
            reference_index: reference_index
            read_group_id: read_group_id
            read_group_fields: read_group_fields
        out:
            [aligned_bam]
    merge:
        run: merge.cwl
        in:
            bams: align/aligned_bam
        out:
            [merged_bam]
