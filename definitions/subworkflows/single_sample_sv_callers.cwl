#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Subworkflow to allow calling different SV callers which require bam files as inputs"

requirements:
    - class: ScatterFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement

inputs:
    bam:
        type: File
        secondaryFiles: [.bai,^.bai]
    reference:
        type: string

    cnvkit_diagram:
        type: boolean?
    cnvkit_drop_low_coverage:
        type: boolean?
    cnvkit_method:
        type: string?
    cnvkit_reference_cnn:
        type: File
    cnvkit_scatter_plot:
        type: boolean?
    cnvkit_male_reference:
        type: boolean?
    cnvkit_vcf_name:
        type: string?

    manta_call_regions:
        type: File?
        secondaryFiles: [.tbi]
    manta_non_wgs:
        type: boolean?
    manta_output_contigs:
        type: boolean?

    merge_max_distance:
        type: int
    merge_min_svs:
        type: int
    merge_same_type:
        type: boolean
    merge_same_strand:
        type: boolean
    merge_estimate_sv_distance:
        type: boolean
    merge_min_sv_size:
        type: int

    smoove_exclude_regions:
        type: File?
    snps_vcf:
        type: File?
    genome_build:
        type: string
    sv_paired_percentage:
        type: double?
    sv_paired_count:
        type: int?
    sv_split_percentage:
        type: double?
    sv_split_count:
        type: int?
    cnv_deletion_depth:
        type: double?
    cnv_duplication_depth:
        type: double?

outputs:
    cn_diagram:
        type: File?
        outputSource: run_cnvkit/cn_diagram
    cn_scatter_plot:
        type: File?
        outputSource: run_cnvkit/cn_scatter_plot
    tumor_antitarget_coverage:
        type: File
        outputSource: run_cnvkit/tumor_antitarget_coverage
    tumor_target_coverage:
        type: File
        outputSource: run_cnvkit/tumor_target_coverage
    tumor_bin_level_ratios:
        type: File
        outputSource: run_cnvkit/tumor_bin_level_ratios
    tumor_segmented_ratios:
        type: File
        outputSource: run_cnvkit/tumor_segmented_ratios
    cnvkit_vcf:
        type: File
        outputSource: run_cnvkit_raw_index/indexed_vcf
    cnvnator_cn_file:
        type: File
        outputSource: run_cnvnator/cn_file
    cnvnator_root:
        type: File
        outputSource: run_cnvnator/root_file
    cnvnator_vcf:
        type: File
        outputSource: run_cnvnator_raw_index/indexed_vcf
    manta_diploid_variants:
        type: File?
        outputSource: run_manta/diploid_variants
    manta_somatic_variants:
        type: File?
        outputSource: run_manta/somatic_variants
    manta_all_candidates:
        type: File
        outputSource: run_manta/all_candidates
    manta_small_candidates:
        type: File
        outputSource: run_manta/small_candidates
    manta_tumor_only_variants:
        type: File?
        outputSource: run_manta/tumor_only_variants
    smoove_output_variants:
        type: File
        outputSource: run_smoove/output_vcf
    filtered_sv_vcfs:
        type: File[]
        outputSource: filtered_vcf_index/indexed_vcf
    merged_sv_vcf:
        type: File
        outputSource: run_merge/merged_sv_vcf
    merged_annotated_tsv:
        type: File
        outputSource: run_merge/merged_annotated_tsv
steps:
    run_cnvkit:
        run: cnvkit_single_sample.cwl
        in:
            diagram: cnvkit_diagram
            drop_low_coverage: cnvkit_drop_low_coverage
            method: cnvkit_method
            reference_cnn: cnvkit_reference_cnn
            tumor_bam: bam
            scatter_plot: cnvkit_scatter_plot
            male_reference: cnvkit_male_reference
            cnvkit_vcf_name: cnvkit_vcf_name
        out:
            [cn_diagram, cn_scatter_plot, tumor_antitarget_coverage, tumor_target_coverage, tumor_bin_level_ratios, tumor_segmented_ratios, cnvkit_vcf]
    run_cnvkit_raw_bgzip:
        run: ../tools/bgzip.cwl
        in:
            file: run_cnvkit/cnvkit_vcf
        out:
            [bgzipped_file]
    run_cnvkit_raw_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: run_cnvkit_raw_bgzip/bgzipped_file
        out:
            [indexed_vcf]
    run_cnvkit_filter:
        run: ../tools/filter_sv_vcf_depth.cwl
        in:
            input_vcf: run_cnvkit/cnvkit_vcf
            deletion_depth: cnv_deletion_depth
            duplication_depth: cnv_duplication_depth
            output_vcf_name:
                default: "filtered_cnvkit.vcf"
            vcf_source:
                default: "cnvkit"
        out:
            [filtered_sv_vcf]
    run_cnvnator:
        run: ../tools/cnvnator.cwl
        in:
            bam: bam
            reference: reference
            sample_name:
                default: "CNVnator"
        out:
            [vcf, root_file, cn_file]
    run_cnvnator_raw_bgzip:
        run: ../tools/bgzip.cwl
        in:
            file: run_cnvnator/vcf
        out:
            [bgzipped_file]
    run_cnvnator_raw_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: run_cnvnator_raw_bgzip/bgzipped_file
        out:
            [indexed_vcf]
    run_cnvnator_filter:
        run: ../tools/filter_sv_vcf_depth.cwl
        in:
            input_vcf: run_cnvnator/vcf
            deletion_depth: cnv_deletion_depth
            duplication_depth: cnv_duplication_depth
            output_vcf_name:
                default: "filtered_CNVnator.vcf"
            vcf_source:
                default: "cnvnator"
        out:
            [filtered_sv_vcf]
    run_manta:
        run: ../tools/manta_somatic.cwl
        in:
            call_regions: manta_call_regions
            non_wgs: manta_non_wgs
            output_contigs: manta_output_contigs
            reference: reference
            tumor_bam: bam
        out:
            [diploid_variants, somatic_variants, all_candidates, small_candidates, tumor_only_variants]
    run_manta_filter:
        run: ../tools/filter_sv_vcf_read_support.cwl
        in:
            input_vcf: run_manta/tumor_only_variants
            output_vcf_name:
                default: "filtered_manta.vcf"
            paired_count: sv_paired_count
            paired_percentage: sv_paired_percentage
            split_count: sv_split_count
            split_percentage: sv_split_percentage
            vcf_source:
                default: "manta"
        out:
            [filtered_sv_vcf]
    run_smoove:
        run: ../tools/smoove.cwl
        in:
            bams:
                source: [bam]
                linkMerge: merge_flattened
            exclude_regions: smoove_exclude_regions
            reference: reference
        out:
            [output_vcf]
    run_smoove_filter:
        run: ../tools/filter_sv_vcf_read_support.cwl
        in:
            input_vcf: run_smoove/output_vcf
            output_vcf_name:
                default: "filtered_smoove.vcf"
            paired_count: sv_paired_count
            paired_percentage: sv_paired_percentage
            split_count: sv_split_count
            split_percentage: sv_split_percentage
            vcf_source:
                default: "smoove"
        out:
            [filtered_sv_vcf]
    run_merge:
        run: merge_svs.cwl
        in:
            estimate_sv_distance: merge_estimate_sv_distance
            genome_build: genome_build
            max_distance_to_merge: merge_max_distance
            minimum_sv_calls: merge_min_svs
            minimum_sv_size: merge_min_sv_size
            same_strand: merge_same_strand
            same_type: merge_same_type
            snps_vcf: snps_vcf
            sv_vcfs:
                source: [run_cnvkit_filter/filtered_sv_vcf, run_cnvnator_filter/filtered_sv_vcf, run_manta_filter/filtered_sv_vcf, run_smoove_filter/filtered_sv_vcf]
                linkMerge: merge_flattened
        out:
            [merged_sv_vcf, merged_annotated_tsv]
    filtered_vcf_bgzip:
        scatter: [file]
        scatterMethod: dotproduct
        run: ../tools/bgzip.cwl
        in:
            file:
                source: [run_cnvkit_filter/filtered_sv_vcf, run_cnvnator_filter/filtered_sv_vcf, run_manta_filter/filtered_sv_vcf, run_smoove_filter/filtered_sv_vcf]
                linkMerge: merge_flattened
        out: [bgzipped_file]
    filtered_vcf_index:
        scatter: [vcf]
        scatterMethod: dotproduct
        run: ../tools/index_vcf.cwl
        in:
            vcf: filtered_vcf_bgzip/bgzipped_file
        out:
            [indexed_vcf]
