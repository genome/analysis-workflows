#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "alignment"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference_index:
        type: File
        secondaryFiles: [".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"]
    instrument_data_bam:
        type: File
outputs:
    aligned_bam:
        type: File
        outputSource: hisat2_align/aligned_bam
steps:
    bam_to_fastq:
        run: bam_to_fastq.cwl
        in:
            bam: instrument_data_bam
        out:
            [fastq1, fastq2]
    hisat2_align:
        run: hisat2_align.cwl
        in:
            reference_index: reference_index
            fastq1: bam_to_fastq/fastq1
            fastq2: bam_to_fastq/fastq2
        out:
            [aligned_bam]
