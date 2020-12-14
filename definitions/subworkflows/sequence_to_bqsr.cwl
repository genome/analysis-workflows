#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Raw sequence data to BQSR"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
          - $import: ../types/trimming_options.yml
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement

inputs:
    unaligned:
        type: ../types/sequence_data.yml#sequence_data[]
    bqsr_known_sites:
        type: File[]
        secondaryFiles: [.tbi]
        doc: "One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis."
    bqsr_intervals:
        type: string[]?
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
    trimming:
        type:
            - ../types/trimming_options.yml#trimming_options
            - "null"
    final_name:
        type: string
        default: 'final'
outputs:
    final_cram:
        type: File
        outputSource: index_cram/indexed_cram
        secondaryFiles: [.bai, ^.bai]
    mark_duplicates_metrics_file:
        type: File
        outputSource: mark_duplicates_and_sort/metrics_file
steps:
    align:
        scatter: [unaligned]
        scatterMethod: dotproduct
        run: sequence_align_and_tag_adapter.cwl
        in:
            unaligned: unaligned
            reference: reference
            trimming: trimming
        out:
            [aligned_cram]
    merge:
        run: ../tools/merge_crams.cwl
        in:
            crams: align/aligned_cram
            name: final_name
            reference: reference
        out:
            [merged_cram]
    name_sort:
        run: ../tools/name_sort_samtools.cwl
        in:
            cram: merge/merged_cram
            reference: reference
        out:
            [name_sorted_cram]
    mark_duplicates_and_sort:
        run: ../tools/mark_duplicates_and_sort.cwl
        in:
            cram: name_sort/name_sorted_cram
            reference: reference
        out:
            [sorted_cram, metrics_file]
    bqsr:
        run: ../tools/bqsr.cwl
        in:
            reference: reference
            cram: mark_duplicates_and_sort/sorted_cram
            intervals: bqsr_intervals
            known_sites: bqsr_known_sites
        out:
            [bqsr_table]
    apply_bqsr:
        run: ../tools/apply_bqsr.cwl
        in:
            reference: reference
            cram: mark_duplicates_and_sort/sorted_cram
            bqsr_table: bqsr/bqsr_table
            output_name: final_name
        out:
            [bqsr_cram]
    index_cram:
        run: ../tools/index_cram.cwl
        in:
            cram: apply_bqsr/bqsr_cram
        out:
            [indexed_cram]
