#/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Immunotherapy Workflow"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference: string
    reference_index:
        type: File #this requires an extra file with the basename
        secondaryFiles: [".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"]
    reference_annotation:
        type: File
    rna_bams:
        type: File[]
    tumor_bams:
        type: File[]
    normal_bams:
        type: File[]
    rna_readgroups:
        type: string[]
    tumor_readgroups:
        type: string[]
    normal_readgroups:
        type: string[]
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
outputs:
    annotated_vcf:
        type: File
        outputSource: pvacseq/annotated_expression_vcf
    annotated_tsv:
        type: File
        outputSource: pvacseq/annotated_variants_tsv
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
            trimming_max_uncalled: triming_max_uncalled
            trimming_min_readlength: trimming_min_readlength
            kallisto_index: kallisto_index
            gene_transcript_lookup_table: gene_transcript_lookup_table
            strand: strand
            refFlat: refFlat
            ribosomal_intervals: ribosomal_intervals
        out:
            [final_bam, stringtie_transcript_gtf, stringtie_gene_expression_tsv]
    somatic:
        run: somatic_exome.cwl
        in:
            reference: reference
            tumor_bams: tumor_bams
            tumor_readgroups: tumor_readgroups
            tumor_name:
            normal_bams: normal_bams
            normal_readgroups: normal_readgroups
            normal_name:
            mills:
            known_indels:
            dbsnp_vcf:
            bqsr_intervals:
            bait_intervals:
            target_intervals:
            per_base_intervals:
            per_target_intervals:
            summary_intervals:
            omni_vcf:
            picard_metric_accumulation_level:
            qc_minimum_mapping_quality:
            qc_minimum_base_quality:
            interval_list:
            cosmic_vcf:
            panel_of_normals_vcf:
            strelka_cpu_reserved:
            mutect_scatter_count:
            mutect_artifact_detection_mode:
            mutect_max_alt_allele_in_normal_fraction:
            mutect_max_alt_alleles_in_normal_count:
            varscan_strand_filter:
            varscan_min_coverage:
            varscan_min_var_freq:
            varscan_p_value:
            varscan_max_normal_freq:
            pindel_insert_size:
            docm_vcf:
            filter_docm_variants:
            vep_cache_dir:
            synonyms_file:
            annotate_coding_only:
            vep_pick:
            cle_vcf_filter:
            variants_to_table_fields:
            variants_to_table_genotype_fields:
            vep_to_table_fields:
            custom_gnomad_vcf:
            custom_clinvar_vcf:
            manta_call_regions:
            manta_non_wgs:
            manta_output_contigs:
        out:
            [tumor_cram, final_vcf]
    germline:
        run: germline_exome_hla_typing.cwl
        in:
            reference: reference
            bams: normal_bams
            readgroups: normal_readgroups
            mills:
            known_indels:
            dbsnp_vcf:
            bqsr_intervals:
            bait_intervals:
            target_intervals:
            per_base_intervals:
            per_target_intervals:
            summary_intervals:
            omni_vcf:
            picard_metric_accumulation_level:
            emit_reference_confidence:
            gvcf_gq_bands:
            intervals:
            vep_cache_dir:
            synonyms_file:
            annotate_coding_only:
            custom_gnomad_vcf:
            qc_minimum_mapping_quality:
            qc_minimum_base_quality:
            custom_clinvar_vcf:
            optitype_name:
        out:
            [final_vcf]
    phase_vcf:
        run: ../subworkflows/phase_vcf.cwl
        in:
            somatic_vcf: somatic/final_vcf
            germline_vcf: germline/final_vcf
            reference: reference
            reference_dict:
            bam: somatic/tumor_cram
        out:
            [phased_vcf]
    pvacseq:
        run: pvacseq.cwl
        in:
            detect_variants_vcf: somatic/final_vcf
            rnaseq_bam: rnaseq/final_bam
            reference_fasta: reference
                type: string
            readcount_minimum_base_quality:
                type: int?
            readcount_minimum_mapping_quality:
                type: int?
            gene_expression_file:
            transcript_expression_file:
            alleles:
                type: string[]
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
            phased_proximal_variants_vcf: phase_vcf/phased_vcf
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
            n_threads:
                type: int?
            variants_to_table_fields:
                type: string[]?
                default: [CHROM,POS,ID,REF,ALT]
            variants_to_table_genotype_fields:
                type: string[]?
                default: [GT,AD,AF,DP,RAD,RAF,RDP,GX,TX]
            vep_to_table_fields:
                type: string[]?
                default: [HGVSc,HGVSp]
        out:
            [annotated_vcf, annotated_tsv, mhc_i_all_epitopes, mhc_i_filtered_epitopes, mhc_i_ranked_epitopes, mhc_ii_all_epitopes, mhc_ii_filtered_epitopes, mhc_ii_ranked_epitopes, combined_all_epitopes, combined_filtered_epitopes, combined_ranked_epitopes]
