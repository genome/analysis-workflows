#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "umi molecular alignment fastq workflow"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
inputs:
    read1_fastq:
        type: File[]
    read2_fastq:
        type: File[]
    sample_name:
        type: string
    library_name:
        type: string[]
    platform_unit:
        type: string[]
    platform:
        type: string[]
    read_structure:
        type: string[]
    reference:
        type: string
    target_intervals:
       type: File?
outputs:
    aligned_cram:
        type: File
        secondaryFiles: [.crai, ^.crai]
        outputSource: alignment_workflow/aligned_cram
    adapter_histogram:
        type: File[]
        outputSource: alignment_workflow/adapter_histogram
    duplex_seq_metrics:
        type: File[]
        outputSource: alignment_workflow/duplex_seq_metrics
steps:
    fastq_to_bam:
        scatter: [read1_fastq, read2_fastq, library_name, platform_unit, platform]
        scatterMethod: dotproduct
        run: ../tools/fastq_to_bam.cwl
        in:
            read1_fastq: read1_fastq
            read2_fastq: read2_fastq
            sample_name: sample_name
            library_name: library_name
            platform_unit: platform_unit
            platform: platform
        out:
            [bam]
    alignment_workflow:
        run: ../subworkflows/molecular_alignment.cwl
        in:
            bam: fastq_to_bam/bam
            sample_name: sample_name
            read_structure: read_structure
            reference: reference
            target_intervals: target_intervals
        out:
            [aligned_cram, adapter_histogram, duplex_seq_metrics]
