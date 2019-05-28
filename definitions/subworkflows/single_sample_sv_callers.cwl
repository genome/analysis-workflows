#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Subworkflow to allow calling different SV callers which require bam files as inputs"

requirements:
    - class: ScatterFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml

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
    merge_sv_pop_freq_db:
        type: File

    smoove_exclude_regions:
        type: File?

    maximum_sv_pop_freq:
        type: float?
        default: 0.1
    sv_filter_interval_lists:
        type: ../types/labelled_file.yml#labelled_file[]

    vep_cache_dir:
        type: string
    vep_ensembl_assembly:
        type: string
        doc: "genome assembly to use in vep. Examples: GRCh38 or GRCm38"
    vep_ensembl_version:
        type: string
        doc: "ensembl version - Must be present in the cache directory. Example: 95"
    vep_ensembl_species:
        type: string
        doc: "ensembl species - Must be present in the cache directory. Examples: homo_sapiens or mus_musculus"
    variants_to_table_fields:
        type: string[]?
        default: [CHROM,POS,ID,REF,ALT,SVLEN,CHR2,END,POPFREQ_AF,POPFREQ_VarID,NSAMP]
    variants_to_table_genotype_fields:
        type: string[]?
        default: [GT]
    vep_to_table_fields:
        type: string[]?
        default: [SYMBOL]
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
        outputSource: run_cnvkit/cnvkit_vcf
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
    merged_annotated_svs:
        type: File
        outputSource: run_merge/merged_annotated_vcf
    sv_pop_filtered_vcf:
        type: File
        outputSource: filter_vcf/sv_pop_filtered_vcf
    filtered_vcfs:
        type: File[]
        outputSource: annotated_filter_index/indexed_vcf
    annotated_tsvs:
        type: File[]
        outputSource: add_vep_fields_to_table/annotated_variants_tsv
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
    run_merge:
        run: merge_svs.cwl
        in:
            vcfs:
                source: [run_cnvkit/cnvkit_vcf, run_manta/tumor_only_variants, run_smoove/output_vcf]
                linkMerge: merge_flattened
            max_distance_to_merge: merge_max_distance
            minimum_sv_calls: merge_min_svs
            same_type: merge_same_type
            same_strand: merge_same_strand
            estimate_sv_distance: merge_estimate_sv_distance
            minimum_sv_size: merge_min_sv_size
            sv_db: merge_sv_pop_freq_db
            vep_cache_dir: vep_cache_dir
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            reference: reference
        out:
            [merged_annotated_vcf, vep_summary]
    filter_vcf:
        run: filter_sv_vcf.cwl
        in:
            maximum_sv_population_frequency: maximum_sv_pop_freq
            sv_intervals: sv_filter_interval_lists
            vcf: run_merge/merged_annotated_vcf
        out:
            [filtered_vcf, sv_pop_filtered_vcf]
    annotated_filter_bgzip:
        run: ../tools/bgzip.cwl
        scatter: [file]
        scatterMethod: dotproduct
        in:
            file: filter_vcf/filtered_vcf
        out:
            [bgzipped_file]
    annotated_filter_index:
        run: ../tools/index_vcf.cwl
        scatter: [vcf]
        scatterMethod: dotproduct
        in:
            vcf: annotated_filter_bgzip/bgzipped_file
        out:
            [indexed_vcf]
    variants_to_table:
        run: ../tools/variants_to_table.cwl
        scatter: [vcf]
        scatterMethod: dotproduct
        in:
            reference: reference
            vcf: annotated_filter_index/indexed_vcf
            fields: variants_to_table_fields
            genotype_fields: variants_to_table_genotype_fields
        out:
            [variants_tsv]
    add_vep_fields_to_table:
        run: ../tools/add_vep_fields_to_table.cwl
        scatter: [vcf, tsv, prefix]
        scatterMethod: dotproduct
        in:
            vcf: annotated_filter_index/indexed_vcf
            vep_fields: vep_to_table_fields
            tsv: variants_to_table/variants_tsv
            prefix:
                source: sv_filter_interval_lists
                valueFrom: $(self.label)
        out: [annotated_variants_tsv]
