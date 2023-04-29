#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "umi molecular alignment fastq workflow"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
inputs:
    sequence:
        type: ../types/sequence_data.yml#sequence_data[]
        label: "sequence: sequencing data and readgroup information"
        doc: |
          sequence represents the sequencing data as either FASTQs or BAMs with accompanying
          readgroup information. Note that in the @RG field ID and SM are required for FASTQs.
          For BAMs, this pipeline assumes that the RG information is already in the header.
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
    umi_aligned_cram:
        type: File
        secondaryFiles: [.crai, ^.crai]
        outputSource: umi_index_cram/indexed_cram
    grouped_aligned_cram:
        type: File
        secondaryFiles: [.crai, ^.crai]
        outputSource: grouped_index_cram/indexed_cram
    consensus_aligned_cram:
        type: File
        secondaryFiles: [.crai, ^.crai]
        outputSource: consensus_index_cram/indexed_cram
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
        run: ../subworkflows/molecular_alignment.cwl
        in:
            bam: sequence_to_bam/bam
            sample_name: sample_name
            read_structure: read_structure
            reference: reference
            target_intervals: target_intervals
        out:
            [umi_aligned_bam, grouped_aligned_bam, consensus_aligned_bam,  adapter_histogram, duplex_seq_metrics]
    umi_bam_to_cram:
        run: ../tools/bam_to_cram.cwl
        in:
            bam: alignment_workflow/umi_aligned_bam
            reference: reference
        out:
            [cram]
    umi_index_cram:
        run: ../tools/index_cram.cwl
        in: 
            cram: umi_bam_to_cram/cram
        out:
            [indexed_cram]
    grouped_bam_to_cram:
        run: ../tools/bam_to_cram.cwl
        in:
            bam: alignment_workflow/grouped_aligned_bam
            reference: reference
        out:
            [cram]
    grouped_index_cram:
        run: ../tools/index_cram.cwl
        in:
            cram: grouped_bam_to_cram/cram
        out:
            [indexed_cram]
    consensus_bam_to_cram:
        run: ../tools/bam_to_cram.cwl
        in:
            bam: alignment_workflow/consensus_aligned_bam
            reference: reference
        out:
            [cram]
    consensus_index_cram:
         run: ../tools/index_cram.cwl
         in:
            cram: consensus_bam_to_cram/cram
         out:
            [indexed_cram]
