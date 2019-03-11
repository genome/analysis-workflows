#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "merge and annotate svs with population allele freq"
requirements:
    - class: MultipleInputFeatureRequirement
inputs:
    vcfs:
        type: File[]
    max_distance_to_merge:
        type: int
    minimum_sv_calls:
        type: int
    same_type:
        type: int
    same_strand:
        type: int
    estimate_sv_distance:
        type: int
    minimum_sv_size:
        type: int
    cohort_name:
        type: string?
    sv_db:
        type: File
outputs:
    merged_annotated_vcf:
        type: File
        outputSource: add_population_frequency/merged_annotated_vcf
steps:
    merge_vcfs:
        run: ../tools/survivor.cwl
        in: 
            vcfs: vcfs
            max_distance_to_merge: max_distance_to_merge
            minimum_sv_calls: minimum_sv_calls
            same_type: same_type
            same_strand: same_strand
            estimate_sv_distance: estimate_sv_distance
            minimum_sv_size: minimum_sv_size
            cohort_name: cohort_name
        out:
            [merged_vcf]
    add_population_frequency:
        run: ../tools/add_sv_population_frequency.cwl
        in:
            vcf: merge_vcfs/merged_vcf
            sv_db: sv_db
            cohort_name: cohort_name
        out:
            [merged_annotated_vcf]
