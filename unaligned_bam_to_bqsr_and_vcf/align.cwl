#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: 'Unaligned to aligned BAM'
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement
inputs:
    bam:
        type: File
    reference:
        type: File
        secondaryFiles: [.fai, .bwt, .sa, .ann, .amb, .pac, .alt]
    readgroup:
        type: string
outputs:
    tagged_bam:
        type: File
        outputSource: align_and_tag/aligned_bam
steps:
    revert_to_fastq:
        run: revert_to_fastq.cwl
        in:
            bam: bam
        out:
            [fastq, second_end_fastq]
    align_and_tag:
        run: align_and_tag.cwl
        in:
        in:
            readgroup: readgroup
            reference: reference
            fastq: revert_to_fastq/fastq
            fastq2: revert_to_fastq/second_end_fastq
        out:
            [aligned_bam]
