#/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Immunotherapy Workflow"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference_index:
        type: File #this requires an extra file with the basename
        secondaryFiles: [".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"]
    reference_annotation:
        type: File
    rna_bams:
        type: File[]
    rna_readgroups:
        type: string[]
    read_group_fields:
        type:
            type: array
            items:
                type: array
                items: string 
    sample_name:
        type: string
    trimming_adapters:
        type: File
    trimming_adapter_trim_end:
        type: string
    trimming_adapter_min_overlap:
        type: int
    trimming_max_uncalled:
        type: int
    trimming_min_readlength:
        type: int
    kallisto_index:
       type: File
    gene_transcript_lookup_table:
       type: File
    strand:
       type: string?
    refFlat:
        type: File
    ribosomal_intervals:
        type: File

    reference: string
    tumor_bams:
        type: File[]
    tumor_readgroups:
        type: string[]
    tumor_name:
        type: string?
        default: 'tumor'
    normal_bams:
        type: File[]
    normal_readgroups:
        type: string[]
    normal_name:
        type: string?
        default: 'normal'
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
        default: false
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
    filter_docm_variants:
        type: boolean?
        default: true
    vep_cache_dir:
        type: string
    synonyms_file:
        type: File?
    annotate_coding_only:
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
    custom_clinvar_vcf:
        type: File?
        secondaryFiles: [.tbi]
    vep_assembly:
        type: string
        doc: Used to explicitly define which assembly to use; required if there are two or more in the same directory
    manta_call_regions:
        type: File?
        secondaryFiles: [.tbi]
    manta_non_wgs:
        type: boolean?
        default: true
    manta_output_contigs:
        type: boolean?
    somalier_vcf:
        type: File

    emit_reference_confidence:
        type: string
    gvcf_gq_bands:
        type: string[]
    gatk_haplotypecaller_intervals:
        type:
            type: array
            items:
                type: array
                items: string
    optitype_name:
        type: string?

    reference_dict:
        type: File?

    readcount_minimum_base_quality:
        type: int?
    readcount_minimum_mapping_quality:
        type: int?
    prediction_algorithms:
        type: string[]
    epitope_lengths:
        type: int[]?
    binding_threshold:
        type: int?
    allele_specific_binding_thresholds:
        type: boolean?
    minimum_fold_change:
        type: float?
    peptide_sequence_length:
        type: int?
    top_score_metric:
        type:
            - "null"
            - type: enum
              symbols: ["lowest", "median"]
    additional_report_columns:
        type:
            - "null"
            - type: enum
              symbols: ["sample_name"]
    fasta_size:
        type: int?
    downstream_sequence_length:
        type: string?
    exclude_nas:
        type: boolean?
    phased_proximal_variants_vcf:
        type: File?
        secondaryFiles: ['.tbi']
    maximum_transcript_support_level:
        type:
            - "null"
            - type: enum
              symbols: ["1", "2", "3", "4", "5"]
    normal_cov:
        type: int?
    tdna_cov:
        type: int?
    trna_cov:
        type: int?
    normal_vaf:
        type: float?
    tdna_vaf:
        type: float?
    trna_vaf:
        type: float?
    expn_val:
        type: float?
    net_chop_method:
        type:
            - "null"
            - type: enum
              symbols: ["cterm", "20s"]
    net_chop_threshold:
        type: float?
    netmhc_stab:
        type: boolean?
    pvacseq_threads:
        type: int?

outputs:
    annotated_vcf:
        type: File
        outputSource: pvacseq/annotated_vcf
    annotated_tsv:
        type: File
        outputSource: pvacseq/annotated_tsv
    mhc_i_all_epitopes:
        type: File?
        outputSource: pvacseq/mhc_i_all_epitopes
    mhc_i_filtered_epitopes:
        type: File?
        outputSource: pvacseq/mhc_i_filtered_epitopes
    mhc_i_ranked_epitopes:
        type: File?
        outputSource: pvacseq/mhc_i_ranked_epitopes
    mhc_ii_all_epitopes:
        type: File?
        outputSource: pvacseq/mhc_ii_all_epitopes
    mhc_ii_filtered_epitopes:
        type: File?
        outputSource: pvacseq/mhc_ii_filtered_epitopes
    mhc_ii_ranked_epitopes:
        type: File?
        outputSource: pvacseq/mhc_ii_ranked_epitopes
    combined_all_epitopes:
        type: File?
        outputSource: pvacseq/combined_all_epitopes
    combined_filtered_epitopes:
        type: File?
        outputSource: pvacseq/combined_filtered_epitopes
    combined_ranked_epitopes:
        type: File?
        outputSource: pvacseq/combined_ranked_epitopes
steps:
    rnaseq:
        run: rnaseq.cwl
        in:
            reference_index: reference_index
            reference_annotation: reference_annotation
            instrument_data_bams: rna_bams
            read_group_id: rna_readgroups
            read_group_fields: read_group_fields
            sample_name: sample_name
            trimming_adapters: trimming_adapters
            trimming_adapter_trim_end: trimming_adapter_trim_end
            trimming_adapter_min_overlap: trimming_adapter_min_overlap
            trimming_max_uncalled: trimming_max_uncalled
            trimming_min_readlength: trimming_min_readlength
            kallisto_index: kallisto_index
            gene_transcript_lookup_table: gene_transcript_lookup_table
            strand: strand
            refFlat: refFlat
            ribosomal_intervals: ribosomal_intervals
        out:
            [final_bam, stringtie_transcript_gtf, stringtie_gene_expression_tsv, transcript_abundance_tsv, transcript_abundance_h5, gene_abundance, metrics, chart, fusion_evidence]
    somatic:
        run: somatic_exome.cwl
        in:
            reference: reference
            tumor_bams: tumor_bams
            tumor_readgroups: tumor_readgroups
            tumor_name: tumor_name
            normal_bams: normal_bams
            normal_readgroups: normal_readgroups
            normal_name: normal_name
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
            filter_docm_variants: filter_docm_variants
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
            vep_assembly: vep_assembly
            manta_call_regions: manta_call_regions
            manta_non_wgs: manta_non_wgs
            manta_output_contigs: manta_output_contigs
            somalier_vcf: somalier_vcf
        out:
            [tumor_cram,tumor_mark_duplicates_metrics,tumor_insert_size_metrics,tumor_alignment_summary_metrics,tumor_hs_metrics,tumor_per_target_coverage_metrics,tumor_per_target_hs_metrics,tumor_per_base_coverage_metrics,tumor_per_base_hs_metrics,tumor_summary_hs_metrics,tumor_flagstats,tumor_verify_bam_id_metrics,tumor_verify_bam_id_depth,normal_cram,normal_mark_duplicates_metrics,normal_insert_size_metrics,normal_alignment_summary_metrics,normal_hs_metrics,normal_per_target_coverage_metrics,normal_per_target_hs_metrics,normal_per_base_coverage_metrics,normal_per_base_hs_metrics,normal_summary_hs_metrics,normal_flagstats,normal_verify_bam_id_metrics,normal_verify_bam_id_depth,mutect_unfiltered_vcf,mutect_filtered_vcf,strelka_unfiltered_vcf,strelka_filtered_vcf,varscan_unfiltered_vcf,varscan_filtered_vcf,pindel_unfiltered_vcf,pindel_filtered_vcf,docm_filtered_vcf,final_vcf,final_filtered_vcf,final_tsv,vep_summary,tumor_snv_bam_readcount_tsv,tumor_indel_bam_readcount_tsv,normal_snv_bam_readcount_tsv,normal_indel_bam_readcount_tsv,intervals_antitarget,intervals_target,normal_antitarget_coverage,normal_target_coverage,reference_coverage,cn_diagram,cn_scatter_plot,tumor_antitarget_coverage,tumor_target_coverage,tumor_bin_level_ratios,tumor_segmented_ratios,diploid_variants,somatic_variants,all_candidates,small_candidates,tumor_only_variants,somalier_concordance_metrics,somalier_concordance_statistics]
    germline:
        run: germline_exome_hla_typing.cwl
        in:
            reference: reference
            bams: normal_bams
            readgroups: normal_readgroups
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
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: gatk_haplotypecaller_intervals
            vep_cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            custom_gnomad_vcf: custom_gnomad_vcf
            qc_minimum_mapping_quality: qc_minimum_mapping_quality
            qc_minimum_base_quality: qc_minimum_base_quality
            custom_clinvar_vcf: custom_clinvar_vcf
            optitype_name: optitype_name
            vep_assembly: vep_assembly
        out:
            [cram,mark_duplicates_metrics,insert_size_metrics,insert_size_histogram,alignment_summary_metrics,hs_metrics,per_target_coverage_metrics,per_target_hs_metrics,per_base_coverage_metrics,per_base_hs_metrics,summary_hs_metrics,flagstats,verify_bam_id_metrics,verify_bam_id_depth,gvcf,final_vcf,coding_vcf,limited_vcf,vep_summary,optitype_tsv,optitype_plot]
    phase_vcf:
        run: ../subworkflows/phase_vcf.cwl
        in:
            somatic_vcf: somatic/final_vcf
            germline_vcf: germline/final_vcf
            reference: reference
            reference_dict: reference_dict
            bam: somatic/tumor_cram
        out:
            [phased_vcf]
    pvacseq:
        run: ../subworkflows/pvacseq.cwl
        in:
            detect_variants_vcf: somatic/final_vcf
            sample_name: sample_name #should this be the same as the above sample_name
            normal_sample_name: normal_name #is this right?
            rnaseq_bam: rnaseq/final_bam
            reference_fasta: reference
            readcount_minimum_base_quality: readcount_minimum_base_quality #should this be different from qc_minimum_base_quality?
            readcount_minimum_mapping_quality: readcount_minimum_mapping_quality #should this be different from qc_minimum_mapping_quality
            gene_expression_file: rnaseq/stringtie_gene_expression_tsv #not sure if this connection is correct
            transcript_expression_file: rnaseq/transcript_abundance_tsv #not sure if this connection is correct
            allele_file: germline/optitype_tsv
            prediction_algorithms: prediction_algorithms
            epitope_lengths: epitope_lengths
            binding_threshold: binding_threshold
            allele_specific_binding_thresholds: allele_specific_binding_thresholds
            minimum_fold_change: minimum_fold_change
            peptide_sequence_length: peptide_sequence_length
            top_score_metric: top_score_metric
            additional_report_columns: additional_report_columns
            fasta_size: fasta_size
            downstream_sequence_length: downstream_sequence_length
            exclude_nas: exclude_nas
            phased_proximal_variants_vcf: phase_vcf/phased_vcf
            maximum_transcript_support_level: maximum_transcript_support_level
            normal_cov: normal_cov
            tdna_cov: tdna_cov
            trna_cov: trna_cov
            normal_vaf: normal_vaf
            tdna_vaf: tdna_vaf
            trna_vaf: trna_vaf
            expn_val: expn_val
            net_chop_method: net_chop_method
            net_chop_threshold: net_chop_threshold
            netmhc_stab: netmhc_stab
            n_threads: pvacseq_threads
            variants_to_table_fields: variants_to_table_fields #the defaults set in this file will override the defaults for this parameter set in the pvacseq subworkflow; is this ok?
            variants_to_table_genotype_fields: variants_to_table_genotype_fields #the defaults set in this file will override the defaults for this parameter set in the pvacseq subworkflow; is this ok?
            vep_to_table_fields: vep_to_table_fields
        out:
            [annotated_vcf, annotated_tsv, mhc_i_all_epitopes, mhc_i_filtered_epitopes, mhc_i_ranked_epitopes, mhc_ii_all_epitopes, mhc_ii_filtered_epitopes, mhc_ii_ranked_epitopes, combined_all_epitopes, combined_filtered_epitopes, combined_ranked_epitopes]
