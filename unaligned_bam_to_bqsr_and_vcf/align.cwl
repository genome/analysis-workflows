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
        secondaryFiles: [.fai, .bwt, .sa, .ann, .amb, .pac]
    readgroup:
        type: string
outputs:
    tagged_sam:
        type: File
        outputSource: tag/tagged_sam
steps:
    revert_to_fastq:
        run: revert_to_fastq.cwl
        in:
            bam: bam
        out:
            [fastq, second_end_fastq]
    bwa_mem:
        run: bwa_mem.cwl
        in:
            readgroup: readgroup
            reference: reference
            fastq: revert_to_fastq/fastq
            fastq2: revert_to_fastq/second_end_fastq
        out:
            [aligned_bam]
    tag:
        run: tag.cwl
        in:
            sam: bwa_mem/aligned_bam
        out:
            [tagged_sam]

