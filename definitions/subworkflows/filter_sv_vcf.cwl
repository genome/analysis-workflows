#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "SV filtering workflow"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
    - class: StepInputExpressionRequirement
inputs:
    maximum_sv_population_frequency:
        type: float
    sv_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    vcf:
        type: File
outputs:
    filtered_vcf:
        type: File[]
        outputSource: intersect_variants/intersect_result
steps:
    filter_sv_pop:
        run: ../tools/filter_vcf_max_sv_pop_freq.cwl
        in:
            vcf: vcf
            maximum_sv_population_frequency: maximum_sv_population_frequency
        out:
            [filtered_vcf]
    interval_to_bed:
        run: ../tools/intervals_to_bed.cwl
        scatter: [interval_list]
        scatterMethod: dotproduct
        in:
            interval_list:
                source: sv_intervals
                valueFrom: $(self.file)
        out:    
            [interval_bed]
    intersect_variants:
        run: ../tools/bedtools_intersect.cwl
        scatter: [file_b, output_name]
        scatterMethod: dotproduct
        in:
            file_a: filter_sv_pop/filtered_vcf
            file_b:
                source: interval_to_bed/interval_bed
            output_name:
                source: sv_intervals
                valueFrom: "intersect-$(self.label).vcf"
        out:    
            [intersect_result]
