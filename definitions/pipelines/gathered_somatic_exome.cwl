#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "gathered exome alignment and somatic variant detection"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    reference: string
    tumor_bams:
        type: File[]
    tumor_readgroups:
        type: string[]
    tumor_cram_name:
        type: string?
        default: 'tumor.cram'
    normal_bams:
        type: File[]
    normal_readgroups:
        type: string[]
    normal_cram_name:
        type: string?
        default: 'normal.cram'
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
    cosmic_vcf:
        type: File?
        secondaryFiles: [.tbi]
    panel_of_normals_vcf:
        type: File?
        secondaryFiles: [.tbi]
    strelka_cpu_reserved:
        type: int?
        default: 8
    mutect_scatter_count:
        type: int
    mutect_artifact_detection_mode:
        type: boolean
    mutect_max_alt_allele_in_normal_fraction:
        type: float?
    mutect_max_alt_alleles_in_normal_count:
        type: int?
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
    pindel_insert_size:
        type: int
        default: 400
    docm_vcf:
        type: File
        secondaryFiles: [.tbi]
    vep_cache_dir:
        type: string?
    synonyms_file:
        type: File?
    annotate_coding_only:
        type: boolean?
    hgvs_annotation:
        type: boolean?
    vep_pick:
        type:
            - "null"
            - type: enum
              symbols: ["pick", "flag_pick", "pick_allele", "per_gene", "pick_allele_gene", "flag_pick_allele", "flag_pick_allele_gene"]
    cle_vcf_filter:
        type: boolean
        default: false
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
    output_dir: 
        type: string
outputs:
    final_outputs:
        type: string[]
        outputSource: gatherer/gathered_files
steps:
    somatic_exome:
        run: somatic_exome.cwl
        in:
            reference: reference
            tumor_bams: tumor_bams
            tumor_readgroups: tumor_readgroups
            tumor_cram_name: tumor_cram_name
            normal_bams: normal_bams
            normal_readgroups: normal_readgroups
            normal_cram_name: normal_cram_name
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
            qc_minimum_mapping_quality: qc_minimum_mapping_quality
            qc_minimum_base_quality: qc_minimum_base_quality
            interval_list: interval_list
            cosmic_vcf: cosmic_vcf
            panel_of_normals_vcf: panel_of_normals_vcf
            strelka_cpu_reserved: strelka_cpu_reserved
            mutect_scatter_count: mutect_scatter_count
            mutect_artifact_detection_mode: mutect_artifact_detection_mode
            mutect_max_alt_allele_in_normal_fraction: mutect_max_alt_allele_in_normal_fraction
            mutect_max_alt_alleles_in_normal_count: mutect_max_alt_alleles_in_normal_count
            varscan_strand_filter: varscan_strand_filter
            varscan_min_coverage: varscan_min_coverage
            varscan_min_var_freq: varscan_min_var_freq
            varscan_p_value: varscan_p_value
            varscan_max_normal_freq: varscan_max_normal_freq
            pindel_insert_size: pindel_insert_size
            docm_vcf: docm_vcf
            vep_cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            hgvs_annotation: hgvs_annotation
            vep_pick: vep_pick
            cle_vcf_filter: cle_vcf_filter
            variants_to_table_fields: variants_to_table_fields
            variants_to_table_genotype_fields: variants_to_table_genotype_fields
            vep_to_table_fields: vep_to_table_fields
            custom_gnomad_vcf: custom_gnomad_vcf
        out:
            [tumor_cram, tumor_mark_duplicates_metrics, tumor_insert_size_metrics, tumor_alignment_summary_metrics, tumor_hs_metrics, tumor_per_target_coverage_metrics, tumor_per_base_coverage_metrics, tumor_per_base_hs_metrics, tumor_summary_hs_metrics, tumor_flagstats, tumor_verify_bam_id_metrics, tumor_verify_bam_id_depth, normal_cram, normal_mark_duplicates_metrics, normal_insert_size_metrics, normal_alignment_summary_metrics, normal_hs_metrics, normal_per_target_coverage_metrics, normal_per_target_hs_metrics, normal_per_base_coverage_metrics, normal_per_base_hs_metrics, normal_summary_hs_metrics, normal_flagstats, normal_verify_bam_id_metrics, normal_verify_bam_id_depth, mutect_unfiltered_vcf, mutect_filtered_vcf, strelka_unfiltered_vcf, strelka_filtered_vcf, varscan_unfiltered_vcf, varscan_filtered_vcf, pindel_unfiltered_vcf, pindel_filtered_vcf, docm_unfiltered_vcf, docm_filtered_vcf, final_vcf, final_filtered_vcf, final_tsv, vep_summary, tumor_snv_bam_readcount_tsv, tumor_indel_bam_readcount_tsv, normal_snv_bam_readcount_tsv, normal_indel_bam_readcount_tsv]
    gatherer:
        run: ../tools/gatherer.cwl
        in:
            output_dir: output_dir
            all_files:
                source: [somatic_exome/tumor_cram, somatic_exome/tumor_mark_duplicates_metrics, somatic_exome/tumor_insert_size_metrics, somatic_exome/tumor_alignment_summary_metrics, somatic_exome/tumor_hs_metrics, somatic_exome/tumor_per_target_coverage_metrics, somatic_exome/tumor_per_base_coverage_metrics, somatic_exome/tumor_per_base_hs_metrics, somatic_exome/tumor_summary_hs_metrics, somatic_exome/tumor_flagstats, somatic_exome/tumor_verify_bam_id_metrics, somatic_exome/tumor_verify_bam_id_depth, somatic_exome/normal_cram, somatic_exome/normal_mark_duplicates_metrics, somatic_exome/normal_insert_size_metrics, somatic_exome/normal_alignment_summary_metrics, somatic_exome/normal_hs_metrics, somatic_exome/normal_per_target_coverage_metrics, somatic_exome/normal_per_target_hs_metrics, somatic_exome/normal_per_base_coverage_metrics, somatic_exome/normal_per_base_hs_metrics, somatic_exome/normal_summary_hs_metrics, somatic_exome/normal_flagstats, somatic_exome/normal_verify_bam_id_metrics, somatic_exome/normal_verify_bam_id_depth, somatic_exome/mutect_unfiltered_vcf, somatic_exome/mutect_filtered_vcf, somatic_exome/strelka_unfiltered_vcf, somatic_exome/strelka_filtered_vcf, somatic_exome/varscan_unfiltered_vcf, somatic_exome/varscan_filtered_vcf, somatic_exome/pindel_unfiltered_vcf, somatic_exome/pindel_filtered_vcf, somatic_exome/docm_unfiltered_vcf, somatic_exome/docm_filtered_vcf, somatic_exome/final_vcf, somatic_exome/final_filtered_vcf, somatic_exome/final_tsv, somatic_exome/vep_summary, somatic_exome/tumor_snv_bam_readcount_tsv, somatic_exome/tumor_indel_bam_readcount_tsv, somatic_exome/normal_snv_bam_readcount_tsv, somatic_exome/normal_indel_bam_readcount_tsv]
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
