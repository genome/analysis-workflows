#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "take bisulfite sequence data through trimming, alignment, and markdup"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
          - $import: ../types/trimming_options.yml
inputs:
    sequence:
        type: ../types/sequence_data.yml#sequence_data[]
        doc: "the unaligned sequence data with readgroup information"
    trimming_options:
        type:
            - ../types/trimming_options.yml#trimming_options
            - "null"
    reference_index:
        type: string
    sample_name:
        type: string
outputs:
    aligned_bam:
        type: File
        outputSource: index_bam/indexed_bam
steps:
    trim_and_align:
        scatter: [sequence]
        scatterMethod: dotproduct
        run: sequence_to_bisulfite_alignment_adapter.cwl
        in:
            sequence: sequence
            trimming_options: trimming_options
            reference_index: reference_index
        out:
            [aligned_bam]
    merge:
        run: ../tools/merge_bams_samtools.cwl
        in:
            bams: trim_and_align/aligned_bam
            name: sample_name
        out:
            [merged_bam]
    name_sort:
        run: ../tools/name_sort.cwl
        in:
            bam: merge/merged_bam
        out:
            [name_sorted_bam]
    biscuit_markdup:
        run: ../tools/biscuit_markdup.cwl
        in:
           bam: name_sort/name_sorted_bam
        out:
            [markdup_bam]
    index_bam:
        run: ../tools/index_bam.cwl
        in:
            bam: biscuit_markdup/markdup_bam
        out:
            [indexed_bam]
