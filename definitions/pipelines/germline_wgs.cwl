#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "wgs alignment and germline variant detection"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
          - $import: ../types/sequence_data.yml
          - $import: ../types/trimming_options.yml
          - $import: ../types/vep_custom_annotation.yml
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
    sequence:
        type: ../types/sequence_data.yml#sequence_data[]
        label: "sequence: sequencing data and readgroup information"
        doc: |
          sequence represents the sequencing data as either FASTQs or BAMs with accompanying
          readgroup information. Note that in the @RG field ID and SM are required.
    trimming:
        type:
            - ../types/trimming_options.yml#trimming_options
            - "null"
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
    picard_metric_accumulation_level:
        type: string
    gvcf_gq_bands:
        type: string[]
    intervals:
        type:
            type: array
            items:
                type: array
                items: string
    ploidy:
        type: int?
    qc_intervals:
        type: File
    variant_reporting_intervals:
        type: File
    vep_cache_dir:
        type:
            - string
            - Directory
    vep_ensembl_assembly:
        type: string
        doc: "genome assembly to use in vep. Examples: GRCh38 or GRCm38"
    vep_ensembl_version:
        type: string
        doc: "ensembl version - Must be present in the cache directory. Example: 95"
    vep_ensembl_species:
        type: string
        doc: "ensembl species - Must be present in the cache directory. Examples: homo_sapiens or mus_musculus"
    vep_plugins:
        type: string[]?
        doc: "array of plugins to use when running vep"
    synonyms_file:
        type: File?
    annotate_coding_only:
        type: boolean?
    bqsr_known_sites:
        type: File[]
        secondaryFiles: [.tbi]
        doc: "One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis."
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
    vep_custom_annotations:
        type: ../types/vep_custom_annotation.yml#vep_custom_annotation[]
        doc: "custom type, check types directory for input format"
    cnvkit_diagram:
        type: boolean?
    cnvkit_drop_low_coverage: 
        type: boolean?
    cnvkit_method:
        type:
          - "null"
          - type: enum
            symbols: ["hybrid", "amplicon", "wgs"]
        default: "wgs"
    cnvkit_reference_cnn: 
        type: File?
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
    sv_filter_alt_abundance_percentage:
        type: double?
    sv_filter_paired_count:
        type: int?
    sv_filter_split_count:
        type: int?
    cnv_filter_deletion_depth:
        type: double?
    cnv_filter_duplication_depth:
        type: double?
    variants_to_table_fields:
         type: string[]?
    variants_to_table_genotype_fields:
         type: string[]?
    vep_to_table_fields:
         type: string[]?
    cnv_filter_min_size:
         type: int?
    blocklist_bedpe:
        type: File?
    disclaimer_text:
        type: string?
        default: 'Workflow source can be found at https://github.com/genome/analysis-workflows'
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
    raw_vcf:
        type: File
        outputSource: detect_variants/raw_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: index_disclaimer_final_vcf/indexed_vcf
        secondaryFiles: [.tbi]
    filtered_vcf:
        type: File
        outputSource: index_disclaimer_filtered_vcf/indexed_vcf
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
    cn_diagram:
        type: File?
        outputSource: sv_detect_variants/cn_diagram
    cn_scatter_plot:
        type: File?
        outputSource: sv_detect_variants/cn_scatter_plot
    tumor_antitarget_coverage:
        type: File
        outputSource: sv_detect_variants/tumor_antitarget_coverage
    tumor_target_coverage:
        type: File
        outputSource: sv_detect_variants/tumor_target_coverage
    tumor_bin_level_ratios:
        type: File
        outputSource: sv_detect_variants/tumor_bin_level_ratios
    tumor_segmented_ratios:
        type: File
        outputSource: sv_detect_variants/tumor_segmented_ratios
    cnvkit_vcf:
        type: File
        outputSource: sv_detect_variants/cnvkit_vcf
    cnvnator_cn_file:
        type: File
        outputSource: sv_detect_variants/cnvnator_cn_file
    cnvnator_root:
        type: File
        outputSource: sv_detect_variants/cnvnator_root
    cnvnator_vcf:
        type: File
        outputSource: sv_detect_variants/cnvnator_vcf
    manta_diploid_variants:
        type: File?
        outputSource: sv_detect_variants/manta_diploid_variants
    manta_somatic_variants:
        type: File?
        outputSource: sv_detect_variants/manta_somatic_variants
    manta_all_candidates:
        type: File
        outputSource: sv_detect_variants/manta_all_candidates
    manta_small_candidates:
        type: File
        outputSource: sv_detect_variants/manta_small_candidates
    manta_tumor_only_variants:
        type: File?
        outputSource: sv_detect_variants/manta_tumor_only_variants
    smoove_output_variants:
        type: File
        outputSource: sv_detect_variants/smoove_output_variants
    final_tsv:
        type: File
        outputSource: add_disclaimer_final_tsv/output_file
    filtered_tsv:
        type: File
        outputSource: add_disclaimer_filtered_tsv/output_file
    cnvkit_filtered_vcf:
        type: File
        outputSource: sv_detect_variants/cnvkit_filtered_vcf
    cnvnator_filtered_vcf:
        type: File
        outputSource: sv_detect_variants/cnvnator_filtered_vcf
    manta_filtered_vcf:
        type: File
        outputSource: sv_detect_variants/manta_filtered_vcf
    smoove_filtered_vcf:
        type: File
        outputSource: sv_detect_variants/smoove_filtered_vcf
    survivor_merged_vcf:
        type: File
        outputSource: add_disclaimer_survivor_sv_vcf/output_file
    survivor_merged_annotated_tsv:
        type: File
        outputSource: add_disclaimer_survivor_sv_tsv/output_file
    bcftools_merged_vcf:
        type: File
        outputSource: add_disclaimer_bcftools_sv_vcf/output_file
    bcftools_merged_annotated_tsv:
        type: File
        outputSource: add_disclaimer_bcftools_sv_tsv/output_file
    bcftools_merged_filtered_annotated_tsv:
        type: File
        outputSource: add_disclaimer_bcftools_filtered_sv_tsv/output_file
steps:
    alignment_and_qc:
        run: alignment_wgs.cwl
        in:
            reference: reference
            sequence: sequence
            trimming: trimming
            omni_vcf: omni_vcf
            intervals: qc_intervals
            picard_metric_accumulation_level: picard_metric_accumulation_level
            bqsr_known_sites: bqsr_known_sites
            bqsr_intervals: bqsr_intervals
            minimum_mapping_quality: minimum_mapping_quality
            minimum_base_quality: minimum_base_quality
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
        out:
            [bam, mark_duplicates_metrics, insert_size_metrics, insert_size_histogram, alignment_summary_metrics, gc_bias_metrics, gc_bias_metrics_chart, gc_bias_metrics_summary, wgs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth, per_base_coverage_metrics, per_base_hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, summary_hs_metrics]
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
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            ploidy: ploidy
            contamination_fraction: extract_freemix/freemix_score
            vep_cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            limit_variant_intervals: variant_reporting_intervals
            vep_custom_annotations: vep_custom_annotations
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            vep_plugins: vep_plugins
            vep_to_table_fields: vep_to_table_fields
            variants_to_table_fields: variants_to_table_fields
            variants_to_table_genotype_fields: variants_to_table_genotype_fields
        out:
            [raw_vcf, final_vcf, filtered_vcf, vep_summary, final_tsv, filtered_tsv]
    add_disclaimer_filtered_vcf:
        run: ../tools/add_string_at_line_bgzipped.cwl
        in:
            input_file: detect_variants/filtered_vcf
            line_number:
                default: 2
            some_text:
                source: disclaimer_text
                valueFrom: "##disclaimer=$(self)"
            output_name:
                source: detect_variants/filtered_vcf
                valueFrom: "$(self.basename)"
        out:
            [output_file]
    index_disclaimer_filtered_vcf:
        run: ../tools/index_vcf.cwl
        in:
            vcf: add_disclaimer_filtered_vcf/output_file
        out:
            [indexed_vcf]
    add_disclaimer_final_vcf:
        run: ../tools/add_string_at_line_bgzipped.cwl
        in:
            input_file: detect_variants/final_vcf
            line_number:
                default: 2
            some_text:
                source: disclaimer_text
                valueFrom: "##disclaimer=$(self)"
            output_name:
                source: detect_variants/final_vcf
                valueFrom: "$(self.basename)"
        out:
            [output_file]
    index_disclaimer_final_vcf:
        run: ../tools/index_vcf.cwl
        in:
            vcf: add_disclaimer_final_vcf/output_file
        out:
            [indexed_vcf]
    add_disclaimer_filtered_tsv:
        run: ../tools/add_string_at_line.cwl
        in:
            input_file: detect_variants/filtered_tsv
            line_number:
                default: 1
            some_text:
                source: disclaimer_text
                valueFrom: "#$(self)"
            output_name:
                source: detect_variants/filtered_tsv
                valueFrom: "$(self.basename)"
        out:
            [output_file]
    add_disclaimer_final_tsv:
        run: ../tools/add_string_at_line.cwl
        in:
            input_file: detect_variants/final_tsv
            line_number:
                default: 1
            some_text:
                source: disclaimer_text
                valueFrom: "#$(self)"
            output_name:
                source: detect_variants/final_tsv
                valueFrom: "$(self.basename)"
        out:
            [output_file]
    sv_detect_variants:
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
            cnv_deletion_depth: cnv_filter_deletion_depth
            cnv_duplication_depth: cnv_filter_duplication_depth
            cnv_filter_min_size: cnv_filter_min_size
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
            snps_vcf: detect_variants/final_vcf
            sv_alt_abundance_percentage: sv_filter_alt_abundance_percentage
            sv_paired_count: sv_filter_paired_count
            sv_split_count: sv_filter_split_count
            genome_build: vep_ensembl_assembly
            blocklist_bedpe: blocklist_bedpe
        out: 
           [cn_diagram, cn_scatter_plot, tumor_antitarget_coverage, tumor_target_coverage, tumor_bin_level_ratios, tumor_segmented_ratios, cnvkit_vcf, cnvnator_cn_file, cnvnator_root, cnvnator_vcf, manta_diploid_variants, manta_somatic_variants, manta_all_candidates, manta_small_candidates, manta_tumor_only_variants, smoove_output_variants, cnvkit_filtered_vcf, cnvnator_filtered_vcf, manta_filtered_vcf, smoove_filtered_vcf, survivor_merged_vcf, survivor_merged_annotated_tsv, bcftools_merged_vcf, bcftools_merged_annotated_tsv, bcftools_merged_filtered_annotated_tsv]
    add_disclaimer_survivor_sv_vcf:
        run: ../tools/add_string_at_line_bgzipped.cwl
        in:
            input_file: sv_detect_variants/survivor_merged_vcf
            line_number:
                default: 2
            some_text:
                source: disclaimer_text
                valueFrom: "##disclaimer=$(self)"
            output_name:
                source: sv_detect_variants/survivor_merged_vcf
                valueFrom: "$(self.basename)"
        out:
            [output_file]
    add_disclaimer_bcftools_sv_vcf:
        run: ../tools/add_string_at_line_bgzipped.cwl
        in:
            input_file: sv_detect_variants/bcftools_merged_vcf
            line_number:
                default: 2
            some_text:
                source: disclaimer_text
                valueFrom: "##disclaimer=$(self)"
            output_name:
                source: sv_detect_variants/bcftools_merged_vcf
                valueFrom: "$(self.basename)"
        out:
            [output_file]
    add_disclaimer_survivor_sv_tsv:
        run: ../tools/add_string_at_line.cwl
        in:
            input_file: sv_detect_variants/survivor_merged_annotated_tsv
            line_number:
                default: 1
            some_text:
                source: disclaimer_text
                valueFrom: "#$(self)"
            output_name:
                source: sv_detect_variants/survivor_merged_annotated_tsv
                valueFrom: "$(self.basename)"
        out:
            [output_file]
    add_disclaimer_bcftools_sv_tsv:
        run: ../tools/add_string_at_line.cwl
        in:
            input_file: sv_detect_variants/bcftools_merged_annotated_tsv
            line_number:
                default: 1
            some_text:
                source: disclaimer_text
                valueFrom: "#$(self)"
            output_name:
                source: sv_detect_variants/bcftools_merged_annotated_tsv
                valueFrom: "$(self.basename)"
        out:
            [output_file]
    add_disclaimer_bcftools_filtered_sv_tsv:
        run: ../tools/add_string_at_line.cwl
        in:
            input_file: sv_detect_variants/bcftools_merged_filtered_annotated_tsv
            line_number:
                default: 1
            some_text:
                source: disclaimer_text
                valueFrom: "#$(self)"
            output_name:
                source: sv_detect_variants/bcftools_merged_filtered_annotated_tsv
                valueFrom: "$(self.basename)"
        out:
            [output_file]
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
