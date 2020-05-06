#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "given multiple somatic vcfs from the same individual, merges the vcfs, adds readcounts, and creates a table"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
    - class: MultipleInputFeatureRequirement
inputs:
    vcfs:
        type: [File]
        secondaryFiles: [.tbi]
    tumor_bams:
        type: [File]
        secondaryFiles: [.bai, .crai]
        doc: array of tumor bams or crams
    tumor_sample_names:
        type: [string]
        doc: tumor sample names - assumes same order as the tumor bams
    normal_bam:
        type: File
        secondaryFiles: [.bai, .crai]
        doc: shared normal bam file
    normal_sample_name:
        type: string
        doc: shared normal sample name
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
    merge_vcfs:
        run: merge_somatic_vcfs.cwl
        in:
            vcfs: vcfs
            sample_names: tumor_sample_names
            normal_sample_name: normal_sample_name
        out:
            [merged_vcf]
    add_readcounts:
        run: bam_readcount_multisample.cwl
        in:
            vcf: merge_vcfs/merged_vcf
            bams:
                source: [normal_bam, bams]
                linkMerge: merge_flattened
            sample_name: 
                source: [normal_sample_name, tumor_sample_names]
                linkMerge: merge_flattened
            reference_fasta: reference_fasta
            min_base_quality: min_base_quality
            min_mapping_quality: min_mapping_quality
            data_type: data_type
        out:
            [readcount_vcfs]
