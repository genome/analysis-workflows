#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "bam to trimmed fastqs and HISAT alignments"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: InlineJavascriptRequirement
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
    read_group_fields:
        type:
            type: array
            items: string
    reference_index:
        type: File
        secondaryFiles: [".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"]
    strand:
        type:
          - "null"
          - type: enum
            symbols: ["first", "second", "unstranded"]
    paired_end:
        type:
            type: enum
            symbols: ["true", "false"]
        default: "true"
        doc: 'whether the sequence data is paired-end (for single-end override to false)'
outputs:
    fastqs:
        type: File[]
        outputSource: bam_to_trimmed_fastq/trimmed_fastqs
    aligned_bam:
        type: File
        outputSource: hisat2_align/aligned_bam
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
    hisat2_align:
        run: ../tools/hisat2_align.cwl
        in:
            reference_index: reference_index
            paired_end: paired_end
            fastqs: bam_to_trimmed_fastq/trimmed_fastqs
            read_group_id: read_group_id
            read_group_fields: read_group_fields
            strand: strand
        out:
            [aligned_bam]
