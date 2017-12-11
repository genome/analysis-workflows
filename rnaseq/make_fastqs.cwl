#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "make_fastqs"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    instrument_data_bam:
        type: File
    trimming_adapters:
        type: File
    trimming_adapter_trim_end:
        type: string
    trimming_adapter_min_overlap:
        type: int
    trimming_max_uncalled:
        type: int
    trimming_min_readlength:
        type: int
outputs:
    trimmed_fastqs:
        type: File
        outputSource: [trimmed_fastq1, trimmed_fastq2]
steps:
    bam_to_fastq:
        run: bam_to_fastq.cwl
        in:
            bam: instrument_data_bam
        out:
            [fastq1, fastq2]
    trim_fastq:
        run: trim_fastq.cwl
        in:
            reads1: bam_to_fastq/fastq1
            reads2: bam_to_fastq/fastq2
            adapters: trimming_adapters
            adapter_trim_end: trimming_adapter_trim_end
            adapter_min_overlap: trimming_adapter_min_overlap
            max_uncalled: trimming_max_uncalled
            min_readlength: trimming_min_readlength
        out:
            [trimmed_fastq1, trimmed_fastq2]
