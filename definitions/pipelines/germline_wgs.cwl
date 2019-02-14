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
    synonyms_file:
        type: File?
    coding_only:
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
    manta_non_wgs:
        type: boolean?
    manta_output_contigs:
        type: boolean?
    smoove_exclude_regions:
        type: File?
outputs:
    cram:
        type: File
        outputSource: alignment_and_qc/cram
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
            [cram, mark_duplicates_metrics, insert_size_metrics, insert_size_histogram, alignment_summary_metrics, gc_bias_metrics, gc_bias_metrics_chart, gc_bias_metrics_summary, wgs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth, per_base_coverage_metrics, per_base_hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, summary_hs_metrics, bamcoverage_bigwig]
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
            cram: alignment_and_qc/cram
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            contamination_fraction: extract_freemix/freemix_score
            vep_cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            coding_only: coding_only
            custom_gnomad_vcf: custom_gnomad_vcf
            limit_variant_intervals: variant_reporting_intervals
            custom_clinvar_vcf: custom_clinvar_vcf
        out:
            [gvcf, final_vcf, coding_vcf, limited_vcf, vep_summary]
    variant_callers:
        run: ../subworkflows/single_sample_sv_callers.cwl
        in:
            cram: alignment_and_qc/cram
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
        out: 
           [cn_diagram, cn_scatter_plot, tumor_antitarget_coverage, tumor_target_coverage, tumor_bin_level_ratios, tumor_segmented_ratios, cnvkit_vcf, manta_diploid_variants, manta_somatic_variants, manta_all_candidates, manta_small_candidates, manta_tumor_only_variants, smoove_output_variants] 
