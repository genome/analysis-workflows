#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "bam to trimmed fastqs and biscuit alignments"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    bam:
        type: File
    adapters:
        type: File
    adapter_trim_end:
        type: string
    adapter_min_overlap:
        type: int
    max_uncalled:
        type: int
    min_readlength:
        type: int
    read_group_id:
        type: string
    reference_index:
        type: string
    paired_end:
        type:
            type: enum
            symbols: ["true", "false"]
        default: "true"
        doc: 'whether the sequence data is paired-end (for single-end override to false)'
outputs:
    aligned_bam:
        type: File
        outputSource: biscuit_markdup/markdup_bam
steps:
    bam_to_trimmed_fastq:
        run: bam_to_trimmed_fastq.cwl
        in:
            bam: bam
            paired_end: paired_end
            adapters: adapters
            adapter_trim_end: adapter_trim_end
            adapter_min_overlap: adapter_min_overlap
            max_uncalled: max_uncalled
            min_readlength: min_readlength
        out:
            [trimmed_fastqs]
    biscuit_align:
        run: ../tools/biscuit_align.cwl
        in:
            reference_index: reference_index
            fastqs: bam_to_trimmed_fastq/trimmed_fastqs
            read_group_id: read_group_id
        out:
            [aligned_bam]
    index_bam:
        run: ../tools/index_bam.cwl
        in:
            bam: biscuit_align/aligned_bam
        out:
            [indexed_bam]
    biscuit_markdup:
        run: ../tools/biscuit_markdup.cwl
        in:
           bam: index_bam/indexed_bam
        out:
            [markdup_bam]
