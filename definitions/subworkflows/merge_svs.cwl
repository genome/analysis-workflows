#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "merge and annotate svs with population allele freq and vep"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    vcfs:
        type: File[]
    max_distance_to_merge:
        type: int
    minimum_sv_calls:
        type: int
    same_type:
        type: boolean
    same_strand:
        type: boolean
    estimate_sv_distance:
        type: boolean
    minimum_sv_size:
        type: int
    cohort_name:
        type: string?
    sv_db:
        type: File
    vep_cache_dir:
        type: string
    coding_only:
        type: boolean?
    custom_gnomad_vcf:
        type: File?
    custom_clinvar_vcf:
        type: File?
    reference:
        type: string
    synonyms_file:
        type: File?
    vep_plugins:
        type: string[]?
        default: []
    vep_assembly:
        type: string
outputs:
    merged_annotated_vcf:
        type: File
        outputSource: sort_vcf/sorted_vcf
    vep_summary:
        type: File
        outputSource: annotate_variants/vep_summary
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
    annotate_variants:
        run: ../tools/vep.cwl
        in:
            assembly: vep_assembly
            vcf: add_population_frequency/merged_annotated_vcf
            cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            coding_only: coding_only
            reference: reference
            custom_gnomad_vcf: custom_gnomad_vcf
            custom_clinvar_vcf: custom_clinvar_vcf
            plugins: vep_plugins
            everything:
                default: false
            pick:
                default: "per_gene"
        out:
            [annotated_vcf, vep_summary]
    sort_vcf:
        run: ../tools/sort_vcf.cwl
        in:
            vcf: annotate_variants/annotated_vcf
        out:
            [sorted_vcf]
