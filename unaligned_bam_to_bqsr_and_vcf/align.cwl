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
    tagged_bam:
        type: File
        outputSource: sam_to_bam/bam
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
            [aligned_sam]
    tag:
        run: tag.cwl
        in:
            sam: bwa_mem/aligned_sam
        out:
            [tagged_sam]

