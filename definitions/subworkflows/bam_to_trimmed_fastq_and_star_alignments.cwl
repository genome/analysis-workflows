#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "bam to trimmed fastqs and STAR alignments"
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
    genomeDir:
        type: Directory
    runMode:
        type: string
    outSAMattrRGline:
        type:
            type: array
            items: string
    outFileNamePrefix:
        type: string
    outSAMstrandField:
        type: string
outputs:
    fastqs:
        type: File[]
        outputSource: trim_fastq/fastqs
    aligned_bam:
        type: File
        outputSource: star_align/aligned_bam
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
            [fastqs]
    star_align:
        run: ../tools/star_align.cwl
        in:
            genomeDir: genomeDir
            runMode: runMode
            outSAMattrRGline: outSAMattrRGline
            outSAMstrandField: outSAMstrandField
            outFileNamePrefix: outFileNamePrefix
            fastq1: 
                source: trim_fastq/fastqs
                valueFrom: $(self[0])
            fastq2: 
                source: trim_fastq/fastqs
                valueFrom: $(self[1])
        out:
            [aligned_bam]
