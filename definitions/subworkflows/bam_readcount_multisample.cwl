#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "given a VCF, gets readcounts from each of a list of samples' bam/cram files and adds them to the VCF"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
inputs:
    vcf:
        type: File
        secondaryFiles: [.tbi]
    sample_names:
        type: [string]
    bams:
        type: [File]
        secondaryFiles: [.bai, .crai]
    reference_fasta:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    min_base_quality:
        type: int?
    min_mapping_quality:
        type: int?
    data_type:
        type:
            - type: enum
              symbols: ["DNA", "RNA"]
        doc: for now, this only accepts either "DNA" or "RNA" and assumes it applies to all samples/bams to avoid having to pass in an array
outputs:
    merged_readcount_vcf:
        type: File
        outputSource: ##bam_readcount/indel_bam_readcount_tsv
steps:
    get_readcounts:
        scatter: [bam, sample_name]
        scatterMethod: dotproduct
        run: bam_readcount_to_vcf.cwl
        in:
            vcf: vcf
            bam: bams
            sample_name: sample_names
            reference_fasta: reference_fasta
            min_base_quality: min_base_quality
            min_mapping_quality: min_mapping_quality
            data_type: data_type
        out:
            [readcount_vcfs, sample_names]

    subset_vcfs_by_sample:
        scatter: [vcf, sample_name]
        scatterMethod: dotproduct
        run: tools/subset_vcf_by_sample.cwl
        in:
            vcf: get_readcounts/readcount_vcfs
            sample_names: get_readcounts/sample_names
        out:
            [subset_vcfs]

    merge_vcfs:
        run: ../tools/merge.cwl
        in:
            vcfs: subset_vcfs_by_sample/subset_vcfs
        out:
            [merged_vcf]
