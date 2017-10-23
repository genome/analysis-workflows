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
    read_group_id:
        type: string
    read_group_fields:
        type: string[]
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
    hisat2_align:
        run: hisat2_align.cwl
        in:
            reference_index: reference_index
            fastq1: trim_fastq/trimmed_fastq1
            fastq2: trim_fastq/trimmed_fastq2
            read_group_id: read_group_id
            read_group_fields: read_group_fields
        out:
            [aligned_bam]
