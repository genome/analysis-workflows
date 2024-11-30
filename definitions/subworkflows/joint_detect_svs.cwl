#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "joint detect svs"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    bams:
        type: File[]
        secondaryFiles: [^.bai]
    sample_names:
        type: string[]
    cohort_name:
        type: string
    exclude_regions:
        type: File?
    manta_call_regions:
        type: File?
    manta_output_contigs:
        type: boolean?
    cnvnator_bin_size:
        type: int?
    cnvkit_method:
        type:
          - "null"
          - type: enum
            symbols: ["hybrid", "amplicon", "wgs"]
    cnvkit_reference_cnn:
        type: File?
    cnvkit_segment_filter:
        type:
          - "null"
          - type: enum
            symbols: ["ampdel", "ci", "cn", "sem"]
    filter_del_depth:
        type: double?
    filter_dup_depth:
        type: double?
    filter_paired_count:
        type: int?
    filter_split_count:
        type: int?
    filter_alt_abundance_percentage:
        type: double?
    filter_depth_caller_min_size:
        type: int?
    survivor_estimate_sv_distance:
        type: boolean
    genome_build:
        type: string
    survivor_max_distance_to_merge:
        type: int
    survivor_minimum_sv_calls:
        type: int
    survivor_minimum_sv_size:
        type: int
    survivor_same_strand:
        type: boolean
    survivor_same_type:
        type: boolean
    snps_vcf:
        type: File?
    filter_blocklist_bedpe:
        type: File?
    annotsv_filter_pop_af:
        type: double?
    annotsv_filter_no_CDS:
        type: boolean?
    annotsv_annotations:
        type:
            - string
            - Directory
outputs:
    all_staged:
        type: Directory
        outputSource: stage_all/gathered_directory
steps:
# stage 1, variant calling
    smoove:
        run: ../tools/smoove.cwl
        in:
            bams: bams
            cohort_name: cohort_name
            reference: reference
            exclude_regions: exclude_regions
        out:
            [output_vcf]
    index_smoove:
        run: ../tools/index_vcf.cwl
        in:
            vcf: smoove/output_vcf
        out:
            [indexed_vcf]
    stage_raw_smoove:
        run: ../tools/gather_to_sub_directory.cwl
        in:
            outdir:
                valueFrom: "smoove"
            files:
                source: [index_smoove/indexed_vcf]
                linkMerge: merge_flattened
        out:
            [gathered_directory]
    manta:
        run: ../tools/manta_germline.cwl
        in:
            bams: bams
            reference: reference
            call_regions: manta_call_regions
            output_contigs: manta_output_contigs
        out:
            [diploid_variants, all_candidates, small_candidates, stats]
    stage_raw_manta:
        run: ../tools/gather_to_sub_directory.cwl
        in:
            outdir:
                valueFrom: "manta"
            files:
                source: [manta/diploid_variants, manta/all_candidates, manta/small_candidates]
                linkMerge: merge_flattened
            directory: manta/stats
        out:
            [gathered_directory]
    cnvnator:
        run: joint_cnvnator.cwl
        in:
            reference: reference
            sample_names: sample_names
            bams: bams
            bin_size: cnvnator_bin_size
        out:
            [vcfs, root_files, cn_files]
    stage_raw_cnvnator:
        run: ../tools/gather_to_sub_directory.cwl
        in:
            outdir:
                valueFrom: "cnvnator"
            files:
                source: [cnvnator/vcfs, cnvnator/root_files, cnvnator/cn_files]
                linkMerge: merge_flattened
        out:
            [gathered_directory]
    cnvkit:
        run: joint_cnvkit.cwl
        in:
            sample_names: sample_names
            bams: bams
            reference_fasta: reference
            reference_cnn: cnvkit_reference_cnn
            method: cnvkit_method
            segment_filter: cnvkit_segment_filter
        out:
            [vcfs, cnr, cns]
    stage_raw_cnvkit:
        run: ../tools/gather_to_sub_directory.cwl
        in:
            outdir:
                valueFrom: "cnvkit"
            files:
                source: [cnvkit/vcfs, cnvkit/cnr, cnvkit/cns]
                linkMerge: merge_flattened
        out:
            [gathered_directory]
    stage_raw:
        run: ../tools/gather_to_sub_directory_dirs.cwl
        in:
             outdir:
                 valueFrom: "raw"
             directories:
                 source: [stage_raw_smoove/gathered_directory, stage_raw_manta/gathered_directory, stage_raw_cnvnator/gathered_directory, stage_raw_cnvkit/gathered_directory]
                 linkMerge: merge_flattened
        out:
            [gathered_directory]
# stage 2, filtering
    filter_smoove:
        run: sv_joint_read_caller_filter.cwl
        in:
            reference: reference
            sample_names: sample_names
            bams: bams
            filter_del_depth: filter_del_depth
            filter_dup_depth: filter_dup_depth
            filter_paired_count: filter_paired_count
            filter_split_count: filter_split_count
            filter_alt_abundance_percentage: filter_alt_abundance_percentage
            sv_vcf: index_smoove/indexed_vcf
            vcf_source:
                default: "smoove"
        out:
            [vcfs]
    filter_manta:
        run: sv_joint_read_caller_filter.cwl
        in:
            reference: reference
            sample_names: sample_names
            bams: bams
            filter_del_depth: filter_del_depth
            filter_dup_depth: filter_dup_depth
            filter_paired_count: filter_paired_count
            filter_split_count: filter_split_count
            filter_alt_abundance_percentage: filter_alt_abundance_percentage
            sv_vcf: manta/diploid_variants
            vcf_source:
                default: "manta"
        out:
            [vcfs]
    filter_cnvnator:
        run: sv_joint_depth_caller_filter.cwl
        in:
            reference: reference
            sample_names: sample_names
            bams: bams
            filter_del_depth: filter_del_depth
            filter_dup_depth: filter_dup_depth
            sv_vcfs: cnvnator/vcfs
            vcf_source:
                default: "cnvnator"
            min_sv_size: filter_depth_caller_min_size
        out:
            [vcfs]
    filter_cnvkit:
        run: sv_joint_depth_caller_filter.cwl
        in:
            reference: reference
            sample_names: sample_names
            bams: bams
            filter_del_depth: filter_del_depth
            filter_dup_depth: filter_dup_depth
            sv_vcfs: cnvkit/vcfs
            vcf_source:
                default: "cnvkit"
            min_sv_size: filter_depth_caller_min_size
        out:
            [vcfs]
    stage_filtered:
        run: ../tools/gather_to_sub_directory.cwl
        in:
            outdir:
                valueFrom: "filtered"
            files:
                source: [filter_smoove/vcfs, filter_manta/vcfs, filter_cnvnator/vcfs, filter_cnvkit/vcfs]
                linkMerge: merge_flattened
        out:
            [gathered_directory]
#  stage3, merge+annotate+filter
    merge_svs:
        run: merge_svs.cwl
        in:
            cohort_name: cohort_name
            estimate_sv_distance: survivor_estimate_sv_distance
            genome_build: genome_build
            max_distance_to_merge: survivor_max_distance_to_merge
            minimum_sv_calls: survivor_minimum_sv_calls
            minimum_sv_size: survivor_minimum_sv_size
            same_strand: survivor_same_strand
            same_type: survivor_same_type
            snps_vcf: snps_vcf
            sv_vcfs:
                source: [filter_smoove/vcfs, filter_manta/vcfs, filter_cnvnator/vcfs, filter_cnvkit/vcfs]
                linkMerge: merge_flattened
            blocklist_bedpe: filter_blocklist_bedpe
            filter_pop_af: annotsv_filter_pop_af
            filter_no_CDS: annotsv_filter_no_CDS
            annotsv_annotations: annotsv_annotations
        out:
            [bcftools_merged_sv_vcf, bcftools_merged_annotated_tsv, bcftools_merged_unannotated_tsv, bcftools_merged_filtered_annotated_tsv, survivor_merged_sv_vcf, survivor_merged_annotated_tsv, survivor_merged_unannotated_tsv, survivor_merged_filtered_annotated_tsv]
    stage_merged:
        run: ../tools/gather_to_sub_directory.cwl
        in:
            outdir:
                valueFrom: "merged"
            files:
                source: [merge_svs/bcftools_merged_sv_vcf, merge_svs/bcftools_merged_annotated_tsv, merge_svs/bcftools_merged_unannotated_tsv, merge_svs/bcftools_merged_filtered_annotated_tsv, merge_svs/survivor_merged_sv_vcf, merge_svs/survivor_merged_annotated_tsv, merge_svs/survivor_merged_unannotated_tsv, merge_svs/survivor_merged_filtered_annotated_tsv]
                linkMerge: merge_flattened
        out:
            [gathered_directory]
    stage_all:
        run: ../tools/gather_to_sub_directory_dirs.cwl
        in:
             outdir:
                 valueFrom: "SV_pipeline"
             directories:
                 source: [stage_raw/gathered_directory, stage_filtered/gathered_directory, stage_merged/gathered_directory]
                 linkMerge: merge_flattened
        out:
            [gathered_directory]
