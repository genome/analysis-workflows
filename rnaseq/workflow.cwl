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
    reference_annotation:
        type: File
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
    sample_name:
        type: string
    trimming_adapters:
        type: File
    trimming_adapter_trim_end:
        type: string
    trimming_adapter_min_overlap:
        type: int
    trimming_max_uncalled:
        type: int
    trimming_min_readlength:
        type: int
outputs:
    aligned_bam:
        type: File
        outputSource: merge/merged_bam
    gtf:
        type: File
        outputSource: stringtie/gtf
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
            trimming_adapters: trimming_adapters
            trimming_adapter_trim_end: trimming_adapter_trim_end
            trimming_adapter_min_overlap: trimming_adapter_min_overlap
            trimming_max_uncalled: trimming_max_uncalled
            trimming_min_readlength: trimming_min_readlength
        out:
            [aligned_bam]
    merge:
        run: merge.cwl
        in:
            bams: align/aligned_bam
        out:
            [merged_bam]
    stringtie:
        run: stringtie.cwl
        in:
            bam: merge/merged_bam
            reference_annotation: reference_annotation
            sample_name: sample_name
        out:
            [gtf]
