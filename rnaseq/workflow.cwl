#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "RnaSeq alignment workflow"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
inputs:
    reference_index:
        type: File #this requires an extra file with the basename
        secondaryFiles: [".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"]
    reference_annotation:
        type: File
    instrument_data_bam:
        type: File[]
    read_group_id:
        type: string[]
    read_group_fields:
        type:
            type: array
            items:
                type: array
                items: string
    sample_name:
        type: string
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
    kallisto_index: 
       type: File
    gene_transcript_lookup_table:
       type: File
outputs:
    aligned_bam:
        type: File
        outputBinding:
            glob: hisat2_align/aligned_bam
    gtf:
        type: File
        outputBinding:
            glob:  stringtie/gtf
    transcript_abundance:
        type: File
        outputBinding:
            glob:  kallisto/transcriptQuant.tsv
    gene_abundance:
        type: File
        outputBinding:
            glob:  kallisto/gene_lengths.tsv
    gene_counts:
        type: File
        outputBinding:
            glob:  kallisto/gene_counts.tsv
    gene_lengths:
        type: File
        outputBinding:
            glob:  kallisto/gene_lengths.tsv
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
    kallisto:
        run: kallisto.cwl
        in:
            kallisto_index: kallisto_index
            fastq1: trim_fastq/trimmed_fastq1
            fastq2: trim_fastq/trimmed_fastq2
        out:
            [expression_transcript_table]
    transcript_to_gene:
        run: transcript_to_gene.cwl
        in:
            transcript_table: kallisto/expression_transcript_table
            gene_transcript_lookup_table: gene_transcript_lookup_table
        out:
            [gene_abundance, gene_counts, gene_lengths] 
    hisat2_align:
        run: hisat2_align.cwl
        in:
            reference_index: reference_index
            trimmed_fastq1: trim_fastq/trimmed_fastq1
            trimmed_fastq2: trim_fastq/trimmed_fastq2
            read_group_id: read_group_id
            read_group_fields: read_group_fields
        out:
            [aligned_bam]
    stringtie:
        run: stringtie.cwl
        in:
            bam: hisat2_align/aligned_bam
            reference_annotation: reference_annotation
            sample_name: sample_name
        out:
            [gtf]
