#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Replace legacy AML Trio Assay"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
          - $import: ../types/sequence_data.yml
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    reference: string
    tumor_sequence:
        type: ../types/sequence_data.yml#sequence_data[]
    tumor_name:
        type: string?
        default: 'tumor'
    normal_sequence:
        type: ../types/sequence_data.yml#sequence_data[]
    normal_name:
        type: string?
        default: 'normal'
    followup_sequence:
        type: ../types/sequence_data.yml#sequence_data[]
    followup_name:
        type: string?
        default: 'followup'
    mills:
        type: File
        secondaryFiles: [.tbi]
    known_indels:
        type: File
        secondaryFiles: [.tbi]
    dbsnp_vcf:
        type: File
        secondaryFiles: [.tbi]
    bqsr_intervals:
        type: string[]
    bait_intervals:
        type: File
    target_intervals:
        type: File
    per_base_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    per_target_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    summary_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
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
    variant_reporting_intervals:
        type: File    
    picard_metric_accumulation_level:
        type: string
    qc_minimum_mapping_quality:
        type: int?
        default: 0
    qc_minimum_base_quality:
        type: int?
        default: 0
    interval_list:
        type: File
    strelka_cpu_reserved:
        type: int?
        default: 8
    mutect_scatter_count:
        type: int
    varscan_strand_filter:
        type: int?
        default: 0
    varscan_min_coverage:
        type: int?
        default: 8
    varscan_min_var_freq:
        type: float?
        default: 0.05
    varscan_p_value:
        type: float?
        default: 0.99
    varscan_max_normal_freq:
        type: float?
    pindel_region_file:
        type: File
    pindel_insert_size:
        type: int
        default: 400
    docm_vcf:
        type: File
        secondaryFiles: [.tbi]
    filter_docm_variants:
        type: boolean?
        default: true
    filter_minimum_depth:
        type: int?
        default: 20
    vep_cache_dir:
        type: string
    synonyms_file:
        type: File?
    annotate_coding_only:
        type: boolean?
    germline_coding_only:
        type: boolean?
    vep_pick:
        type:
            - "null"
            - type: enum
              symbols: ["pick", "flag_pick", "pick_allele", "per_gene", "pick_allele_gene", "flag_pick_allele", "flag_pick_allele_gene"]
    cle_vcf_filter:
        type: boolean
        default: true
    variants_to_table_fields:
        type: string[]
        default: [CHROM,POS,ID,REF,ALT,set,AC,AF]
    variants_to_table_genotype_fields:
        type: string[]
        default: [GT,AD]
    vep_to_table_fields:
        type: string[]
        default: [HGVSc,HGVSp]
    vep_ensembl_assembly:
        type: string
        doc: "genome assembly to use in vep. Examples: GRCh38 or GRCm38"
    vep_ensembl_version:
        type: string
        doc: "ensembl version - Must be present in the cache directory. Example: 95"
    vep_ensembl_species:
        type: string
        doc: "ensembl species - Must be present in the cache directory. Examples: homo_sapiens or mus_musculus"
    somalier_vcf:
        type: File
    custom_gnomad_vcf:
        type: File?
        secondaryFiles: [.tbi]
    custom_clinvar_vcf:
        type: File?
        secondaryFiles: [.tbi]
    germline_tsv_prefix:
        type: string?
        default: 'germline_variants'
    germline_variants_to_table_fields:
        type: string[]
        default: [CHROM,POS,ID,REF,ALT,AC,AF]
    germline_variants_to_table_genotype_fields:
        type: string[]
        default: [GT,AD]
    germline_vep_to_table_fields:
        type: string[]
        default: [HGVSc,HGVSp]
    disclaimer_text:
        type: string?
        default: "#This laboratory developed test (LDT) was developed and its performance characteristics determined by the CLIA Licensed Environment laboratory at the McDonnell Genome Institute at Washington University (MGI-CLE, CLIA #26D2092546, CAP #9047655), Dr. David H. Spencer MD, PhD, FCAP, Medical Director. 4444 Forest Park Avenue, Rm 4127 St. Louis, Missouri 63108 (314) 286-1460 Fax: (314) 286-1810. The MGI-CLE laboratory is regulated under CLIA as certified to perform high-complexity testing. This test has not been cleared or approved by the FDA."
    disclaimer_version:
        type: string
outputs:
    tumor_cram:
        type: File
        outputSource: tumor_index_cram/indexed_cram
    tumor_mark_duplicates_metrics:
        type: File
        outputSource: tumor_alignment_and_qc/mark_duplicates_metrics
    tumor_insert_size_metrics:
        type: File
        outputSource: tumor_alignment_and_qc/insert_size_metrics
    tumor_alignment_summary_metrics:
        type: File
        outputSource: tumor_alignment_and_qc/alignment_summary_metrics
    tumor_hs_metrics:
        type: File
        outputSource: tumor_alignment_and_qc/hs_metrics
    tumor_summary_hs_metrics:
        type: File[]
        outputSource: tumor_alignment_and_qc/summary_hs_metrics
    tumor_flagstats:
        type: File
        outputSource: tumor_alignment_and_qc/flagstats
    tumor_verify_bam_id_metrics:
        type: File
        outputSource: tumor_alignment_and_qc/verify_bam_id_metrics
    tumor_verify_bam_id_depth:
        type: File
        outputSource: tumor_alignment_and_qc/verify_bam_id_depth
    normal_cram:
        type: File
        outputSource: normal_index_cram/indexed_cram
    normal_mark_duplicates_metrics:
        type: File
        outputSource: normal_alignment_and_qc/mark_duplicates_metrics
    normal_insert_size_metrics:
        type: File
        outputSource: normal_alignment_and_qc/insert_size_metrics
    normal_alignment_summary_metrics:
        type: File
        outputSource: normal_alignment_and_qc/alignment_summary_metrics
    normal_hs_metrics:
        type: File
        outputSource: normal_alignment_and_qc/hs_metrics
    normal_summary_hs_metrics:
        type: File[]
        outputSource: normal_alignment_and_qc/summary_hs_metrics
    normal_flagstats:
        type: File
        outputSource: normal_alignment_and_qc/flagstats
    normal_verify_bam_id_metrics:
        type: File
        outputSource: normal_alignment_and_qc/verify_bam_id_metrics
    normal_verify_bam_id_depth:
        type: File
        outputSource: normal_alignment_and_qc/verify_bam_id_depth
    followup_cram:
        type: File
        outputSource: followup_index_cram/indexed_cram
    followup_mark_duplicates_metrics:
        type: File
        outputSource: followup_alignment_and_qc/mark_duplicates_metrics
    followup_insert_size_metrics:
        type: File
        outputSource: followup_alignment_and_qc/insert_size_metrics
    followup_alignment_summary_metrics:
        type: File
        outputSource: followup_alignment_and_qc/alignment_summary_metrics
    followup_hs_metrics:
        type: File
        outputSource: followup_alignment_and_qc/hs_metrics
    followup_summary_hs_metrics:
        type: File[]
        outputSource: followup_alignment_and_qc/summary_hs_metrics
    followup_flagstats:
        type: File
        outputSource: followup_alignment_and_qc/flagstats
    followup_verify_bam_id_metrics:
        type: File
        outputSource: followup_alignment_and_qc/verify_bam_id_metrics
    followup_verify_bam_id_depth:
        type: File
        outputSource: followup_alignment_and_qc/verify_bam_id_depth
    mutect_unfiltered_vcf:
        type: File
        outputSource: tumor_detect_variants/mutect_unfiltered_vcf
        secondaryFiles: [.tbi]
    mutect_filtered_vcf:
        type: File
        outputSource: tumor_detect_variants/mutect_filtered_vcf
        secondaryFiles: [.tbi]
    strelka_unfiltered_vcf:
        type: File
        outputSource: tumor_detect_variants/strelka_unfiltered_vcf
        secondaryFiles: [.tbi]
    strelka_filtered_vcf:
        type: File
        outputSource: tumor_detect_variants/strelka_filtered_vcf
        secondaryFiles: [.tbi]
    varscan_unfiltered_vcf:
        type: File
        outputSource: tumor_detect_variants/varscan_unfiltered_vcf
        secondaryFiles: [.tbi]
    varscan_filtered_vcf:
        type: File
        outputSource: tumor_detect_variants/varscan_filtered_vcf
        secondaryFiles: [.tbi]
    pindel_unfiltered_vcf:
        type: File
        outputSource: tumor_detect_variants/pindel_unfiltered_vcf
        secondaryFiles: [.tbi]
    pindel_filtered_vcf:
        type: File
        outputSource: tumor_detect_variants/pindel_filtered_vcf
        secondaryFiles: [.tbi]
    docm_filtered_vcf:
        type: File
        outputSource: tumor_detect_variants/docm_filtered_vcf
        secondaryFiles: [.tbi]
    pindel_region_vcf:
        type: File
        outputSource: pindel_region/pindel_region_vcf
        secondaryFiles: [.tbi]
    tumor_snv_bam_readcount_tsv:
        type: File
        outputSource: tumor_detect_variants/tumor_snv_bam_readcount_tsv
    tumor_indel_bam_readcount_tsv:
        type: File
        outputSource: tumor_detect_variants/tumor_indel_bam_readcount_tsv
    normal_snv_bam_readcount_tsv:
        type: File
        outputSource: tumor_detect_variants/normal_snv_bam_readcount_tsv
    normal_indel_bam_readcount_tsv:
        type: File
        outputSource: tumor_detect_variants/normal_indel_bam_readcount_tsv
    followup_snv_bam_readcount_tsv:
        type: File
        outputSource: followup_bam_readcount/snv_bam_readcount_tsv
    followup_indel_bam_readcount_tsv:
        type: File
        outputSource: followup_bam_readcount/indel_bam_readcount_tsv
    tumor_final_vcf:
        type: File
        outputSource: tumor_detect_variants/final_vcf
        secondaryFiles: [.tbi]
    tumor_final_filtered_vcf:
        type: File
        outputSource: tumor_detect_variants/final_filtered_vcf
        secondaryFiles: [.tbi]
    tumor_final_tsv:
        type: File
        outputSource: add_disclaimer_version_to_tumor_final_tsv/output_file
    tumor_vep_summary:
        type: File
        outputSource: tumor_detect_variants/vep_summary
    germline_final_vcf:
        type: File
        outputSource: germline_detect_variants/final_vcf
        secondaryFiles: [.tbi]
    germline_coding_vcf:
        type: File
        outputSource: germline_detect_variants/coding_vcf
        secondaryFiles: [.tbi]
    germline_limited_vcf:
        type: File
        outputSource: germline_detect_variants/limited_vcf
        secondaryFiles: [.tbi]
    germline_final_tsv:
        type: File
        outputSource: add_disclaimer_version_to_germline_final_tsv/output_file
    somalier_concordance_metrics:
        type: File
        outputSource: concordance/somalier_pairs
    somalier_concordance_statistics:
        type: File
        outputSource: concordance/somalier_samples
    alignment_stat_report:
        type: File
        outputSource: alignment_stat_report/alignment_stat
    coverage_stat_report:
        type: File
        outputSource: coverage_stat_report/coverage_stat
    full_variant_report:
        type: File
        outputSource: add_disclaimer_version_to_full_variant_report/output_file
steps:
    normal_alignment_and_qc:
        run: alignment_exome.cwl
        in:
            reference: reference
            sequence: normal_sequence
            mills: mills
            known_indels: known_indels
            dbsnp_vcf: dbsnp_vcf
            bqsr_intervals: bqsr_intervals
            bait_intervals: bait_intervals
            target_intervals: target_intervals
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
            omni_vcf: omni_vcf
            picard_metric_accumulation_level: picard_metric_accumulation_level   
            minimum_mapping_quality: qc_minimum_mapping_quality
            minimum_base_quality: qc_minimum_base_quality
            final_name:
                source: normal_name
                valueFrom: "$(self).bam"
        out:
            [bam, mark_duplicates_metrics, insert_size_metrics, alignment_summary_metrics, hs_metrics, summary_hs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
    tumor_alignment_and_qc:
        run: alignment_exome.cwl
        in:
            reference: reference
            sequence: tumor_sequence
            mills: mills
            known_indels: known_indels
            dbsnp_vcf: dbsnp_vcf
            bqsr_intervals: bqsr_intervals
            bait_intervals: bait_intervals
            target_intervals: target_intervals
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
            omni_vcf: omni_vcf
            picard_metric_accumulation_level: picard_metric_accumulation_level   
            minimum_mapping_quality: qc_minimum_mapping_quality
            minimum_base_quality: qc_minimum_base_quality
            final_name:
                source: tumor_name
                valueFrom: "$(self).bam"
        out:
            [bam, mark_duplicates_metrics, insert_size_metrics, alignment_summary_metrics, hs_metrics, summary_hs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
    followup_alignment_and_qc:
        run: alignment_exome.cwl
        in:
            reference: reference
            sequence: followup_sequence
            mills: mills
            known_indels: known_indels
            dbsnp_vcf: dbsnp_vcf
            bqsr_intervals: bqsr_intervals
            bait_intervals: bait_intervals
            target_intervals: target_intervals
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
            omni_vcf: omni_vcf
            picard_metric_accumulation_level: picard_metric_accumulation_level   
            minimum_mapping_quality: qc_minimum_mapping_quality
            minimum_base_quality: qc_minimum_base_quality
            final_name:
                source: followup_name
                valueFrom: "$(self).bam"
        out:
            [bam, mark_duplicates_metrics, insert_size_metrics, alignment_summary_metrics, hs_metrics, summary_hs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
    concordance:
        run: ../tools/concordance.cwl
        in:
            reference: reference
            bam_1: tumor_alignment_and_qc/bam
            bam_2: normal_alignment_and_qc/bam
            bam_3: followup_alignment_and_qc/bam
            vcf: somalier_vcf
        out:
            [somalier_pairs, somalier_samples]
    tumor_detect_variants:
        run: detect_variants.cwl
        in:
            reference: reference
            tumor_bam: tumor_alignment_and_qc/bam
            normal_bam: normal_alignment_and_qc/bam
            interval_list: interval_list
            strelka_exome_mode:
                default: true
            strelka_cpu_reserved: strelka_cpu_reserved
            mutect_scatter_count: mutect_scatter_count
            varscan_strand_filter: varscan_strand_filter
            varscan_min_coverage: varscan_min_coverage
            varscan_min_var_freq: varscan_min_var_freq
            varscan_p_value: varscan_p_value
            varscan_max_normal_freq: varscan_max_normal_freq
            pindel_insert_size: pindel_insert_size
            docm_vcf: docm_vcf
            filter_docm_variants: filter_docm_variants
            filter_minimum_depth: filter_minimum_depth
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            vep_cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            vep_pick: vep_pick
            cle_vcf_filter: cle_vcf_filter
            variants_to_table_fields: variants_to_table_fields
            variants_to_table_genotype_fields: variants_to_table_genotype_fields
            vep_to_table_fields: vep_to_table_fields
            custom_gnomad_vcf: custom_gnomad_vcf
            custom_clinvar_vcf: custom_clinvar_vcf
        out:
            [mutect_unfiltered_vcf, mutect_filtered_vcf, strelka_unfiltered_vcf, strelka_filtered_vcf, varscan_unfiltered_vcf, varscan_filtered_vcf, pindel_unfiltered_vcf, pindel_filtered_vcf, docm_filtered_vcf, final_vcf, final_filtered_vcf, final_tsv, vep_summary, tumor_snv_bam_readcount_tsv, tumor_indel_bam_readcount_tsv, normal_snv_bam_readcount_tsv, normal_indel_bam_readcount_tsv]
    add_disclaimer_to_tumor_final_tsv:
        run: ../tools/add_string_at_line.cwl
        in:
            input_file: tumor_detect_variants/final_tsv
            line_number:
                default: 1
            some_text: disclaimer_text
        out:
            [output_file]
    add_disclaimer_version_to_tumor_final_tsv:
        run: ../tools/add_string_at_line.cwl
        in:
            input_file: add_disclaimer_to_tumor_final_tsv/output_file
            line_number:
                default: 2
            some_text: disclaimer_version
        out:
            [output_file]
    pindel_region:
        run: ../subworkflows/pindel_region.cwl
        in:
            reference: reference
            tumor_bam: tumor_alignment_and_qc/bam
            normal_bam: normal_alignment_and_qc/bam
            region_file: pindel_region_file
            insert_size: pindel_insert_size
        out:
            [pindel_region_vcf]
    followup_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: tumor_detect_variants/final_filtered_vcf
            sample:
                default: 'TUMOR'
            reference_fasta: reference
            bam: followup_alignment_and_qc/bam
            prefix: followup_name
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    extract_freemix:
        in:
            verify_bam_id_metrics: normal_alignment_and_qc/verify_bam_id_metrics
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
    germline_detect_variants:
        run: ../subworkflows/germline_detect_variants.cwl
        in:
            reference: reference
            bam: normal_alignment_and_qc/bam
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            contamination_fraction: extract_freemix/freemix_score
            vep_cache_dir: vep_cache_dir
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            synonyms_file: synonyms_file
            annotate_coding_only: germline_coding_only
            custom_gnomad_vcf: custom_gnomad_vcf
            limit_variant_intervals: variant_reporting_intervals
            custom_clinvar_vcf: custom_clinvar_vcf
            variants_to_table_fields: germline_variants_to_table_fields
            variants_to_table_genotype_fields: germline_variants_to_table_genotype_fields
            vep_to_table_fields: germline_vep_to_table_fields
            final_tsv_prefix: germline_tsv_prefix
        out:
            [final_vcf, coding_vcf, limited_vcf, final_tsv]
    add_disclaimer_to_germline_final_tsv:
        run: ../tools/add_string_at_line.cwl
        in:
            input_file: germline_detect_variants/final_tsv
            line_number:
                default: 1
            some_text: disclaimer_text
        out:
            [output_file]
    add_disclaimer_version_to_germline_final_tsv:
        run: ../tools/add_string_at_line.cwl
        in:
            input_file: add_disclaimer_to_germline_final_tsv/output_file
            line_number:
                default: 2
            some_text: disclaimer_version
        out:
            [output_file]
    alignment_stat_report:
        run: ../tools/cle_aml_trio_report_alignment_stat.cwl
        in:
            normal_alignment_summary_metrics: normal_alignment_and_qc/alignment_summary_metrics
            tumor_alignment_summary_metrics: tumor_alignment_and_qc/alignment_summary_metrics
            followup_alignment_summary_metrics: followup_alignment_and_qc/alignment_summary_metrics
        out:
            [alignment_stat]
    coverage_stat_report:
        run: ../tools/cle_aml_trio_report_coverage_stat.cwl
        in:
            normal_roi_hs_metrics: normal_alignment_and_qc/hs_metrics
            normal_summary_hs_metrics: [normal_alignment_and_qc/summary_hs_metrics]
            tumor_roi_hs_metrics: tumor_alignment_and_qc/hs_metrics
            tumor_summary_hs_metrics: [tumor_alignment_and_qc/summary_hs_metrics]
            followup_roi_hs_metrics: followup_alignment_and_qc/hs_metrics
            followup_summary_hs_metrics: [followup_alignment_and_qc/summary_hs_metrics]
        out:
            [coverage_stat]
    full_variant_report:
        run: ../tools/cle_aml_trio_report_full_variants.cwl
        in: 
            variant_tsv: tumor_detect_variants/final_tsv
            followup_snv_bam_readcount: followup_bam_readcount/snv_bam_readcount_tsv
            followup_indel_bam_readcount: followup_bam_readcount/indel_bam_readcount_tsv
            pindel_region_vcf: pindel_region/pindel_region_vcf
        out:
            [full_variant_report]
    add_disclaimer_to_full_variant_report:
        run: ../tools/add_string_at_line.cwl
        in:
            input_file: full_variant_report/full_variant_report
            line_number:
                default: 1
            some_text: disclaimer_text
        out:
            [output_file]
    add_disclaimer_version_to_full_variant_report:
        run: ../tools/add_string_at_line.cwl
        in:
            input_file: add_disclaimer_to_full_variant_report/output_file
            line_number:
                default: 2
            some_text: disclaimer_version
        out:
            [output_file]
    normal_bam_to_cram:
        run: ../tools/bam_to_cram.cwl
        in:
            bam: normal_alignment_and_qc/bam
            reference: reference
        out:
            [cram]
    normal_index_cram:
         run: ../tools/index_cram.cwl
         in:
            cram: normal_bam_to_cram/cram
         out:
            [indexed_cram]
    tumor_bam_to_cram:
        run: ../tools/bam_to_cram.cwl
        in:
            bam: tumor_alignment_and_qc/bam
            reference: reference
        out:
            [cram]
    tumor_index_cram:
         run: ../tools/index_cram.cwl
         in:
            cram: tumor_bam_to_cram/cram
         out:
            [indexed_cram]
    followup_bam_to_cram:
        run: ../tools/bam_to_cram.cwl
        in:
            bam: followup_alignment_and_qc/bam
            reference: reference
        out:
            [cram]
    followup_index_cram:
         run: ../tools/index_cram.cwl
         in:
            cram: followup_bam_to_cram/cram
         out:
            [indexed_cram]
