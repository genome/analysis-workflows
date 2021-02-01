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
    final_bam:
        type: File
        outputSource: index_bam/indexed_bam
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
            [aligned_bam]
    merge:
        run: ../tools/merge_bams.cwl
        in:
            bams: align/aligned_bam
            name: final_name
        out:
            [merged_bam]
    name_sort:
        run: ../tools/name_sort.cwl
        in:
            bam: merge/merged_bam
        out:
            [name_sorted_bam]
    mark_duplicates_and_sort:
        run: ../tools/mark_duplicates_and_sort.cwl
        in:
            bam: name_sort/name_sorted_bam
        out:
            [sorted_bam, metrics_file]
    bqsr:
        run: ../tools/bqsr.cwl
        in:
            reference: reference
            bam: mark_duplicates_and_sort/sorted_bam
            intervals: bqsr_intervals
            known_sites: bqsr_known_sites
        out:
            [bqsr_table]
    apply_bqsr:
        run: ../tools/apply_bqsr.cwl
        in:
            reference: reference
            bam: mark_duplicates_and_sort/sorted_bam
            bqsr_table: bqsr/bqsr_table
            output_name: final_name
        out:
            [bqsr_bam]
    index_bam:
        run: ../tools/index_bam.cwl
        in:
            bam: apply_bqsr/bqsr_bam
        out:
            [indexed_bam]
