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
    QCannotation:
        type: File
#nonmethylation_sites_bed is a variable to determine if user wants to obtain vcf/bed/bedgraph/bigwig for non - methylation sites. Value has to be yes or no
    assay_non_cpg_sites:
        type: string?
        default: 'no'
outputs:
    cram:
        type: File
        outputSource: index_cram/indexed_cram
        secondaryFiles: [.crai, ^.crai]
    vcf:
        type: File
        outputSource: pileup/vcf
    cpgs:
        type: File[]
        outputSource: vcf2bed/final_bed
    cpg_bigwig:
        type: File[]
        outputSource: bedgraph_to_bigwig/final_bigwigs
    gathered_directory:
        type: Directory
        outputSource: bisulfite_qc/QC_directory
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
    bisulfite_qc:
        run: ../subworkflows/bisulfite_qc.cwl
        in:
            vcf: pileup/vcf
            bam: merge/merged_bam
            reference: reference_index
            QCannotation: QCannotation
        out:
            [QC_directory]
    vcf2bed:
        run: ../tools/bisulfite_vcf2bed.cwl
        in:
            vcf: pileup/vcf
            reference: reference_index
            assay_non_cpg_sites: assay_non_cpg_sites
        out:
            [final_bed,final_bedgraph]
    bedgraph_to_bigwig:
        run: ../tools/bedgraph_to_bigwig.cwl
        in:
            final_bedgraph: vcf2bed/final_bedgraph
            reference_sizes: reference_sizes
        out:
            [final_bigwigs]
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
