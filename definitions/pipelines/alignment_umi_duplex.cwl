#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "umi duplex alignment fastq workflow"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
inputs:
    sequence:
        type: ../types/sequence_data.yml#sequence_data[]
    sample_name:
        type: string
    read_structure:
        type: string[]
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
    target_intervals:
       type: File?
outputs:
    aligned_bam:
        type: File
        secondaryFiles: [^.bai]
        outputSource: alignment_workflow/aligned_bam
    adapter_histogram:
        type: File[]
        outputSource: alignment_workflow/adapter_histogram
    duplex_seq_metrics:
        type: File[]
        outputSource: alignment_workflow/duplex_seq_metrics
steps:
    sequence_to_bam:
        scatter: [sequence]
        scatterMethod: dotproduct
        run: ../tools/sequence_to_bam.cwl
        in:
            sequence: sequence
        out:
            [bam]
    alignment_workflow:
        run: ../subworkflows/duplex_alignment.cwl
        in:
            bam: sequence_to_bam/bam
            sample_name: sample_name
            read_structure: read_structure
            reference: reference
            target_intervals: target_intervals
        out:
            [aligned_bam, adapter_histogram, duplex_seq_metrics]
