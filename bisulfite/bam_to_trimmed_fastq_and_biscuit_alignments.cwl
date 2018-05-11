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
outputs:
    aligned_bam:
        type: File
        outputSource: biscuit_align/aligned_bam
steps:
    bam_to_fastq:
        run: bam_to_fastq.cwl
        in:
            bam: bam
        out:
            [fastq1, fastq2]
    trim_fastq:
        run: trim_fastq.cwl
        in:
            reads1: bam_to_fastq/fastq1
            reads2: bam_to_fastq/fastq2
            adapters: adapters
            adapter_trim_end: adapter_trim_end
            adapter_min_overlap: adapter_min_overlap
            max_uncalled: max_uncalled
            min_readlength: min_readlength
        out:
            [fastq1, fastq2]
    biscuit_align:
        run: biscuit_align.cwl
        in:
            reference_index: reference_index
            fastq1: trim_fastq/fastq1
            fastq2: trim_fastq/fastq2
            read_group_id: read_group_id
        out:
            [aligned_bam]
