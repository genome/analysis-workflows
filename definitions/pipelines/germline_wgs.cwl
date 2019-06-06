#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "wgs alignment and germline variant detection"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
    - class: SubworkflowFeatureRequirement
inputs:
    reference: string
    bams:
        type: File[]
    readgroups:
        type: string[]
    mills:
        type: File
        secondaryFiles: [.tbi]
    known_indels:
        type: File
        secondaryFiles: [.tbi]
    dbsnp_vcf:
        type: File
        secondaryFiles: [.tbi]
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
    picard_metric_accumulation_level:
        type: string
    emit_reference_confidence:
        type: string
    gvcf_gq_bands:
        type: string[]
    intervals:
        type:
            type: array
            items:
                type: array
                items: string
    qc_intervals:
        type: File
    variant_reporting_intervals:
        type: File
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
    synonyms_file:
        type: File?
    annotate_coding_only:
        type: boolean?
    custom_gnomad_vcf:
        type: File?
        secondaryFiles: [.tbi]
    bqsr_intervals:
        type: string[]?
    minimum_mapping_quality:
        type: int?
    minimum_base_quality:
        type: int?
    per_base_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    per_target_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    summary_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    custom_clinvar_vcf:
        type: File?
        secondaryFiles: [.tbi]
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
    smoove_exclude_regions:
        type: File?
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
    variants_to_table_fields:
         type: string[]?
    variants_to_table_genotype_fields:
         type: string[]?
    vep_to_table_fields:
         type: string[]?
    maximum_sv_pop_freq:
        type: float?
    sv_filter_interval_lists:
        type: ../types/labelled_file.yml#labelled_file[]
    sv_variants_to_table_fields:
        type: string[]?
    sv_variants_to_table_genotype_fields:
        type: string[]?
    sv_vep_to_table_fields:
        type: string[]?
outputs:
    cram:
        type: File
        outputSource: index_cram/indexed_cram
    mark_duplicates_metrics:
        type: File
        outputSource: alignment_and_qc/mark_duplicates_metrics
    insert_size_metrics:
        type: File
        outputSource: alignment_and_qc/insert_size_metrics
    insert_size_histogram:
        type: File
        outputSource: alignment_and_qc/insert_size_histogram
    alignment_summary_metrics:
        type: File
        outputSource: alignment_and_qc/alignment_summary_metrics
    gc_bias_metrics:
        type: File
        outputSource: alignment_and_qc/gc_bias_metrics
    gc_bias_metrics_chart:
        type: File
        outputSource: alignment_and_qc/gc_bias_metrics_chart
    gc_bias_metrics_summary:
        type: File
        outputSource: alignment_and_qc/gc_bias_metrics_summary
    wgs_metrics:
        type: File
        outputSource: alignment_and_qc/wgs_metrics
    flagstats:
        type: File
        outputSource: alignment_and_qc/flagstats
    verify_bam_id_metrics:
        type: File
        outputSource: alignment_and_qc/verify_bam_id_metrics
    verify_bam_id_depth:
        type: File
        outputSource: alignment_and_qc/verify_bam_id_depth
    gvcf:
        type: File[]
        outputSource: detect_variants/gvcf
    final_vcf:
        type: File
        outputSource: detect_variants/final_vcf
        secondaryFiles: [.tbi]
    coding_vcf:
        type: File
        outputSource: detect_variants/coding_vcf
        secondaryFiles: [.tbi]
    limited_vcf:
        type: File
        outputSource: detect_variants/limited_vcf
        secondaryFiles: [.tbi]
    vep_summary:
        type: File
        outputSource: detect_variants/vep_summary
    per_base_coverage_metrics:
        type: File[]
        outputSource: alignment_and_qc/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File[]
        outputSource: alignment_and_qc/per_base_hs_metrics
    per_target_coverage_metrics:
        type: File[]
        outputSource: alignment_and_qc/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File[]
        outputSource: alignment_and_qc/per_target_hs_metrics
    summary_hs_metrics:
        type: File[]
        outputSource: alignment_and_qc/summary_hs_metrics
    bamcoverage_bigwig:
        type: File
        outputSource: alignment_and_qc/bamcoverage_bigwig
    cn_diagram:
        type: File?
        outputSource: variant_callers/cn_diagram
    cn_scatter_plot:
        type: File?
        outputSource: variant_callers/cn_scatter_plot
    tumor_antitarget_coverage:
        type: File
        outputSource: variant_callers/tumor_antitarget_coverage
    tumor_target_coverage:
        type: File
        outputSource: variant_callers/tumor_target_coverage
    tumor_bin_level_ratios:
        type: File
        outputSource: variant_callers/tumor_bin_level_ratios
    tumor_segmented_ratios:
        type: File
        outputSource: variant_callers/tumor_segmented_ratios
    cnvkit_vcf:
        type: File
        outputSource: variant_callers/cnvkit_vcf
    manta_diploid_variants:
        type: File?
        outputSource: variant_callers/manta_diploid_variants
    manta_somatic_variants:
        type: File?
        outputSource: variant_callers/manta_somatic_variants
    manta_all_candidates:
        type: File
        outputSource: variant_callers/manta_all_candidates
    manta_small_candidates:
        type: File
        outputSource: variant_callers/manta_small_candidates
    manta_tumor_only_variants:
        type: File?
        outputSource: variant_callers/manta_tumor_only_variants
    smoove_output_variants:
        type: File
        outputSource: variant_callers/smoove_output_variants
    merged_annotated_svs:
        type: File
        outputSource: variant_callers/merged_annotated_svs
    final_tsv:
        type: File
        outputSource: detect_variants/final_tsv
    filtered_sv_pop_vcf:
        type: File
        outputSource: variant_callers/sv_pop_filtered_vcf
    filtered_sv_vcfs:
        type: File[]
        outputSource: variant_callers/filtered_vcfs
    annotated_sv_tsvs:
        type: File[]
        outputSource: variant_callers/annotated_tsvs
steps:
    alignment_and_qc:
        run: wgs_alignment.cwl
        in:
            reference: reference
            bams: bams
            readgroups: readgroups
            mills: mills
            known_indels: known_indels
            dbsnp_vcf: dbsnp_vcf
            omni_vcf: omni_vcf
            intervals: qc_intervals
            picard_metric_accumulation_level: picard_metric_accumulation_level
            bqsr_intervals: bqsr_intervals
            minimum_mapping_quality: minimum_mapping_quality
            minimum_base_quality: minimum_base_quality
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
        out:
            [bam, mark_duplicates_metrics, insert_size_metrics, insert_size_histogram, alignment_summary_metrics, gc_bias_metrics, gc_bias_metrics_chart, gc_bias_metrics_summary, wgs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth, per_base_coverage_metrics, per_base_hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, summary_hs_metrics, bamcoverage_bigwig]
    extract_freemix:
        in:
            verify_bam_id_metrics: alignment_and_qc/verify_bam_id_metrics
        out:
            [freemix_score]
        run:
            class: ExpressionTool
            requirements:
                - class: InlineJavascriptRequirement
            inputs:
                verify_bam_id_metrics:
                    type: File
                    inputBinding:
                        loadContents: true
            outputs:
                freemix_score:
                    type: string?
            expression: |
                        ${
                            var metrics = inputs.verify_bam_id_metrics.contents.split("\n");
                            if ( metrics[0].split("\t")[6] == 'FREEMIX' ) {
                                return {'freemix_score': metrics[1].split("\t")[6]};
                            } else {
                                return {'freemix_score:': null };
                            }
                        }
    detect_variants:
        run: ../subworkflows/germline_detect_variants.cwl
        in:
            reference: reference
            bam: alignment_and_qc/bam
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            contamination_fraction: extract_freemix/freemix_score
            vep_cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            custom_gnomad_vcf: custom_gnomad_vcf
            limit_variant_intervals: variant_reporting_intervals
            custom_clinvar_vcf: custom_clinvar_vcf
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            vep_to_table_fields: vep_to_table_fields
            variants_to_table_fields: variants_to_table_fields
            variants_to_table_genotype_fields: variants_to_table_genotype_fields
        out:
            [gvcf, final_vcf, coding_vcf, limited_vcf, vep_summary, final_tsv]
    variant_callers:
        run: ../subworkflows/single_sample_sv_callers.cwl
        in:
            bam: alignment_and_qc/bam
            reference: reference
            cnvkit_diagram: cnvkit_diagram
            cnvkit_drop_low_coverage: cnvkit_drop_low_coverage
            cnvkit_method: cnvkit_method
            cnvkit_reference_cnn: cnvkit_reference_cnn
            cnvkit_scatter_plot: cnvkit_scatter_plot
            cnvkit_male_reference: cnvkit_male_reference
            cnvkit_vcf_name: cnvkit_vcf_name
            manta_call_regions: manta_call_regions
            manta_non_wgs: manta_non_wgs
            manta_output_contigs: manta_output_contigs
            smoove_exclude_regions: smoove_exclude_regions
            merge_max_distance: merge_max_distance
            merge_min_svs: merge_min_svs
            merge_same_type: merge_same_type
            merge_same_strand: merge_same_strand
            merge_estimate_sv_distance: merge_estimate_sv_distance
            merge_min_sv_size: merge_min_sv_size
            merge_sv_pop_freq_db: merge_sv_pop_freq_db
            sv_filter_interval_lists: sv_filter_interval_lists
            vep_cache_dir: vep_cache_dir
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            maximum_sv_pop_freq: maximum_sv_pop_freq
            variants_to_table_fields: sv_variants_to_table_fields
            variants_to_table_genotype_fields: sv_variants_to_table_genotype_fields
            vep_to_table_fields: sv_vep_to_table_fields
        out: 
           [cn_diagram, cn_scatter_plot, tumor_antitarget_coverage, tumor_target_coverage, tumor_bin_level_ratios, tumor_segmented_ratios, cnvkit_vcf, manta_diploid_variants, manta_somatic_variants, manta_all_candidates, manta_small_candidates, manta_tumor_only_variants, smoove_output_variants, merged_annotated_svs, sv_pop_filtered_vcf, filtered_vcfs, annotated_tsvs]
    bam_to_cram:
        run: ../tools/bam_to_cram.cwl
        in:
            bam: alignment_and_qc/bam
            reference: reference
        out:
            [cram]
    index_cram:
         run: ../tools/index_cram.cwl
         in:
            cram: bam_to_cram/cram
         out:
            [indexed_cram]
