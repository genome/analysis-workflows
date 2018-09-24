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
        type: string
    target_intervals:
       type: File?
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
