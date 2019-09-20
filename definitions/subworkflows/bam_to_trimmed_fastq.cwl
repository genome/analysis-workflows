#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "bam to trimmed fastqs"
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
    paired_end:
        type:
            type: enum
            symbols: ["true", "false"]
        default: "true"
        doc: 'whether the sequence data is paired-end (for single-end override to false)'

outputs:
    trimmed_fastqs:
        type: File[]
        outputSource: trim_fastq/trimmed_fastqs
steps:
    bam_to_fastq:
        run: ../tools/bam_to_fastq.cwl
        in:
            bam: bam
            paired_end: paired_end
        out:
            [fastqs]
    trim_fastq:
        run: ../tools/trim_fastq.cwl
        in:
            fastqs: bam_to_fastq/fastqs
            adapters: adapters
            adapter_trim_end: adapter_trim_end
            adapter_min_overlap: adapter_min_overlap
            max_uncalled: max_uncalled
            min_readlength: min_readlength
        out:
            [trimmed_fastqs]
