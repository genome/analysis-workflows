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
    
outputs:
    fastqs:
        type: File[]
        outputSource: trim_fastq/fastqs
    fastq1:
         type: File
         outputSource: trim_fastq/fastq1
    fastq2:
         type: File
         outputSource: trim_fastq/fastq2

steps:
    bam_to_fastq:
        run: ../tools/bam_to_fastq.cwl
        in:
            bam: bam
        out:
            [fastq1, fastq2]
    trim_fastq:
        run: ../tools/trim_fastq.cwl
        in:
            reads1: bam_to_fastq/fastq1
            reads2: bam_to_fastq/fastq2
            adapters: adapters
            adapter_trim_end: adapter_trim_end
            adapter_min_overlap: adapter_min_overlap
            max_uncalled: max_uncalled
            min_readlength: min_readlength
        out:
            [fastqs, fastq1, fastq2]
    
