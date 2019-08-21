#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "umi duplex alignment workflow"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
inputs:
    bam:
        type: File[]
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
    min_reads:
       type: int[]
       default: [10, 5, 3]
    max_read_error_rate:
       type: float?
       default: 0.05
    max_base_error_rate:
       type: float?
       default: 0.1
    min_base_quality:
       type: int
       default: 50
    max_no_call_fraction:
       type: float
       default: 0.05
outputs:
    aligned_bam:
        type: File
        secondaryFiles: [^.bai]
        outputSource: clip_overlap/clipped_bam
    adapter_histogram:
        type: File[]
        outputSource: align/adapter_metrics
    duplex_seq_metrics:
        type: File[]
        outputSource: collect_duplex_seq_metrics/duplex_seq_metrics
steps:
    align:
        scatter: bam
        run: umi_alignment.cwl
        in:
            bam: bam
            read_structure: read_structure
            reference: reference
        out:
            [aligned_bam, adapter_metrics]
    merge:
        run: ../tools/merge_bams.cwl
        in:
            bams: align/aligned_bam
        out:
            [merged_bam]
    group_reads_by_umi:
        run: ../tools/group_reads.cwl
        in:
            bam: merge/merged_bam
        out:
            [grouped_bam]
    call_duplex_consensus:
        run: ../tools/call_duplex_consensus.cwl
        in:
            bam: group_reads_by_umi/grouped_bam
        out:
            [consensus_bam]
    align_consensus:
        run: ../tools/realign.cwl
        in:
            bam: call_duplex_consensus/consensus_bam
            reference: reference
        out:
            [consensus_aligned_bam]
    filter_consensus:
        run: ../tools/filter_consensus.cwl
        in:
            bam: align_consensus/consensus_aligned_bam
            reference: reference
            min_reads: min_reads
            max_read_error_rate: max_read_error_rate
            max_base_error_rate: max_base_error_rate
            min_base_quality: min_base_quality
            max_no_call_fraction: max_no_call_fraction
        out:
            [filtered_bam]
    clip_overlap:
        run: ../tools/clip_overlap.cwl
        in:
            bam: filter_consensus/filtered_bam
            reference: reference
        out:
            [clipped_bam]
    collect_duplex_seq_metrics:
       run: ../tools/duplex_seq_metrics.cwl
       in:
            bam: group_reads_by_umi/grouped_bam
            intervals: target_intervals
            description: sample_name
       out:
            [duplex_seq_metrics]
