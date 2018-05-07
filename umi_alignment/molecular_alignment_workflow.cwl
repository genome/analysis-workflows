#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "umi alignment workflow"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    bam:
        type: File
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
        type: File
        outputSource: mark_illumina_adapters/metrics
    duplex_seq_metrics:
        type: File[]
        outputSource: collect_duplex_seq_metrics/duplex_seq_metrics
steps:
    extract_umis:
        run: extract_umis.cwl
        in:
            bam: bam
            read_structure: read_structure
        out:
            [umi_extracted_bam]
    mark_illumina_adapters:
        run: mark_illumina_adapters.cwl
        in:
            bam: extract_umis/umi_extracted_bam
        out:
            [marked_bam, metrics]
    align:
        run: align.cwl
        in:
            bam: mark_illumina_adapters/marked_bam
            reference: reference
        out:
            [aligned_bam]
    group_reads_by_umi:
        run: group_reads.cwl
        in:
            bam: align/aligned_bam
        out:
            [grouped_bam]
    call_molecular_consensus:
        run: call_molecular_consensus.cwl
        in:
            bam: group_reads_by_umi/grouped_bam
        out:
            [consensus_bam]
    align_consensus:
        run: realign.cwl
        in:
            bam: call_molecular_consensus/consensus_bam
            reference: reference
        out:
            [consensus_aligned_bam]
    filter_consensus:
        run: filter_consensus.cwl
        in:
            bam: align_consensus/consensus_aligned_bam
            reference: reference
        out:
            [filtered_bam]
    clip_overlap:
        run: clip_overlap.cwl
        in:
            bam: filter_consensus/filtered_bam
            reference: reference
        out:
            [clipped_bam]
    collect_duplex_seq_metrics:
       run: duplex_seq_metrics.cwl
       in:
            bam: group_reads_by_umi/grouped_bam
            intervals: target_intervals
            description: sample_name
       out:
            [duplex_seq_metrics]
