#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "gather AML trio outputs"
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
    somalier_vcf:
        type: File
    output_dir: 
        type: string
    disclaimer_version:
        type: string
outputs:
    final_outputs:
        type: string[]
        outputSource: gatherer/gathered_files
steps:
    aml_trio:
        run: aml_trio_cle.cwl
        in:
            reference: reference
            tumor_sequence: tumor_sequence
            tumor_name: tumor_name
            normal_sequence: normal_sequence
            normal_name: normal_name
            followup_sequence: followup_sequence
            followup_name: followup_name
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
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            variant_reporting_intervals: variant_reporting_intervals
            picard_metric_accumulation_level: picard_metric_accumulation_level   
            qc_minimum_mapping_quality: qc_minimum_mapping_quality
            qc_minimum_base_quality: qc_minimum_base_quality
            interval_list: interval_list
            strelka_cpu_reserved: strelka_cpu_reserved
            mutect_scatter_count: mutect_scatter_count
            varscan_strand_filter: varscan_strand_filter
            varscan_min_coverage: varscan_min_coverage
            varscan_min_var_freq: varscan_min_var_freq
            varscan_p_value: varscan_p_value
            varscan_max_normal_freq: varscan_max_normal_freq
            pindel_region_file: pindel_region_file
            pindel_insert_size: pindel_insert_size
            docm_vcf: docm_vcf
            filter_docm_variants: filter_docm_variants
            filter_minimum_depth: filter_minimum_depth
            vep_cache_dir: vep_cache_dir
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            germline_coding_only: germline_coding_only
            vep_pick: vep_pick
            cle_vcf_filter: cle_vcf_filter
            variants_to_table_fields: variants_to_table_fields
            variants_to_table_genotype_fields: variants_to_table_genotype_fields
            vep_to_table_fields: vep_to_table_fields
            custom_gnomad_vcf: custom_gnomad_vcf
            custom_clinvar_vcf: custom_clinvar_vcf
            somalier_vcf: somalier_vcf
            germline_tsv_prefix: germline_tsv_prefix
            germline_variants_to_table_fields: germline_variants_to_table_fields
            germline_variants_to_table_genotype_fields: germline_variants_to_table_genotype_fields
            germline_vep_to_table_fields: germline_vep_to_table_fields 
            disclaimer_version: disclaimer_version
        out:
            [tumor_cram, tumor_mark_duplicates_metrics, tumor_insert_size_metrics, tumor_alignment_summary_metrics, tumor_hs_metrics, tumor_summary_hs_metrics, tumor_flagstats, tumor_verify_bam_id_metrics, tumor_verify_bam_id_depth, normal_cram, normal_mark_duplicates_metrics, normal_insert_size_metrics, normal_alignment_summary_metrics, normal_hs_metrics, normal_summary_hs_metrics, normal_flagstats, normal_verify_bam_id_metrics, normal_verify_bam_id_depth, followup_cram, followup_mark_duplicates_metrics, followup_insert_size_metrics, followup_alignment_summary_metrics, followup_hs_metrics, followup_summary_hs_metrics, followup_flagstats, followup_verify_bam_id_metrics, followup_verify_bam_id_depth, mutect_unfiltered_vcf, mutect_filtered_vcf, strelka_unfiltered_vcf, strelka_filtered_vcf, varscan_unfiltered_vcf, varscan_filtered_vcf, pindel_unfiltered_vcf, pindel_filtered_vcf, docm_filtered_vcf, pindel_region_vcf, tumor_final_vcf, tumor_final_filtered_vcf, tumor_final_tsv, tumor_vep_summary, germline_final_vcf, germline_coding_vcf, germline_limited_vcf, germline_final_tsv, alignment_stat_report, coverage_stat_report, full_variant_report, tumor_snv_bam_readcount_tsv, tumor_indel_bam_readcount_tsv, normal_snv_bam_readcount_tsv, normal_indel_bam_readcount_tsv, followup_snv_bam_readcount_tsv, followup_indel_bam_readcount_tsv, somalier_concordance_metrics, somalier_concordance_statistics]
    gatherer:
        run: ../tools/gatherer.cwl
        in:
            output_dir: output_dir
            all_files:
                source: [aml_trio/tumor_cram, aml_trio/tumor_mark_duplicates_metrics, aml_trio/tumor_insert_size_metrics, aml_trio/tumor_alignment_summary_metrics, aml_trio/tumor_hs_metrics, aml_trio/tumor_summary_hs_metrics, aml_trio/tumor_flagstats, aml_trio/tumor_verify_bam_id_metrics, aml_trio/tumor_verify_bam_id_depth, aml_trio/normal_cram, aml_trio/normal_mark_duplicates_metrics, aml_trio/normal_insert_size_metrics, aml_trio/normal_alignment_summary_metrics, aml_trio/normal_hs_metrics, aml_trio/normal_summary_hs_metrics, aml_trio/normal_flagstats, aml_trio/normal_verify_bam_id_metrics, aml_trio/normal_verify_bam_id_depth, aml_trio/followup_cram, aml_trio/followup_mark_duplicates_metrics, aml_trio/followup_insert_size_metrics, aml_trio/followup_alignment_summary_metrics, aml_trio/followup_hs_metrics, aml_trio/followup_summary_hs_metrics, aml_trio/followup_flagstats, aml_trio/followup_verify_bam_id_metrics, aml_trio/followup_verify_bam_id_depth, aml_trio/mutect_unfiltered_vcf, aml_trio/mutect_filtered_vcf, aml_trio/strelka_unfiltered_vcf, aml_trio/strelka_filtered_vcf, aml_trio/varscan_unfiltered_vcf, aml_trio/varscan_filtered_vcf, aml_trio/pindel_unfiltered_vcf, aml_trio/pindel_filtered_vcf, aml_trio/docm_filtered_vcf, aml_trio/pindel_region_vcf, aml_trio/tumor_final_vcf, aml_trio/tumor_final_filtered_vcf, aml_trio/tumor_final_tsv, aml_trio/tumor_vep_summary, aml_trio/germline_final_vcf, aml_trio/germline_coding_vcf, aml_trio/germline_limited_vcf, aml_trio/germline_final_tsv, aml_trio/alignment_stat_report, aml_trio/coverage_stat_report, aml_trio/full_variant_report, aml_trio/tumor_snv_bam_readcount_tsv, aml_trio/tumor_indel_bam_readcount_tsv, aml_trio/normal_snv_bam_readcount_tsv, aml_trio/normal_indel_bam_readcount_tsv, aml_trio/followup_snv_bam_readcount_tsv, aml_trio/followup_indel_bam_readcount_tsv, aml_trio/somalier_concordance_metrics, aml_trio/somalier_concordance_statistics]
                valueFrom: ${
                                function flatten(inArr, outArr) {
                                    var arrLen = inArr.length;
                                    for (var i = 0; i < arrLen; i++) {
                                        if (Array.isArray(inArr[i])) {
                                            flatten(inArr[i], outArr);
                                        }
                                        else {
                                            outArr.push(inArr[i]);
                                        }
                                    }
                                    return outArr;
                                }
                                var no_secondaries = flatten(self, []);
                                var all_files = []; 
                                var arrLen = no_secondaries.length;
                                for (var i = 0; i < arrLen; i++) {
                                    all_files.push(no_secondaries[i]);
                                    var secondaryLen = no_secondaries[i].secondaryFiles.length;
                                    for (var j = 0; j < secondaryLen; j++) {
                                        all_files.push(no_secondaries[i].secondaryFiles[j]);
                                    }
                                }
                                return all_files;
                            }
        out: [gathered_files]
