#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Bisulfite alignment and QC"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
inputs:
    reference_index:
        type: string
    reference_sizes:
        type: File
    instrument_data_bams:
        type: File[]
    read_group_id:
        type: string[]
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
outputs:
    cram:
        type: File
        outputSource: index_cram/indexed_cram
        secondaryFiles: [.crai, ^.crai]
    vcf:
        type: File
        outputSource: pileup/vcf
    cpgs:
        type: File
        outputSource: vcf2bed/cpgs
    cpg_bedgraph:
        type: File
        outputSource: bedgraph_to_bigwig/cpg_bigwig
steps:
    bam_to_trimmed_fastq_and_biscuit_alignments:
        run: ../subworkflows/bam_to_trimmed_fastq_and_biscuit_alignments.cwl
        scatter: [bam, read_group_id]
        scatterMethod: dotproduct
        in:
            bam: instrument_data_bams
            read_group_id: read_group_id
            adapters: trimming_adapters
            adapter_trim_end: trimming_adapter_trim_end
            adapter_min_overlap: trimming_adapter_min_overlap
            max_uncalled: trimming_max_uncalled
            min_readlength: trimming_min_readlength
            reference_index: reference_index
        out:
            [aligned_bam]
    merge:
        run: ../tools/merge_bams.cwl
        in:
            bams: bam_to_trimmed_fastq_and_biscuit_alignments/aligned_bam
        out:
            [merged_bam]
    pileup:
        run: ../tools/biscuit_pileup.cwl
        in:
            bam: merge/merged_bam
            reference: reference_index
        out:
            [vcf]
    vcf2bed:
        run: ../tools/bisulfite_vcf2bed.cwl
        in:
            vcf: pileup/vcf
            reference: reference_index
        out:
            [cpgs,cpg_bedgraph]
    bedgraph_to_bigwig:
        run: ../tools/bedgraph_to_bigwig.cwl
        in:
            bedgraph: vcf2bed/cpg_bedgraph
            reference_sizes: reference_sizes
        out:
            [cpg_bigwig]
    bam_to_cram:
        run: ../tools/bam_to_cram.cwl
        in:
            reference: reference_index
            bam: merge/merged_bam
        out:
            [cram]
    index_cram:
        run: ../tools/index_cram.cwl
        in:
            cram: bam_to_cram/cram
        out:
            [indexed_cram]
