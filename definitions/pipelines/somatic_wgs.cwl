#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Whole genome alignment and somatic variant detection"
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
    reference: string
    tumor_sequence:
        type: ../types/sequence_data.yml#sequence_data[]
        label: "tumor_sequence: MT sequencing data and readgroup information"
        doc: |
          tumor_sequence represents the sequencing data for the MT sample as either FASTQs or BAMs with
          accompanying readgroup information. Note that in the @RG field ID and SM are required.
    tumor_name:
        type: string?
        default: 'tumor'
    normal_sequence:
        type: ../types/sequence_data.yml#sequence_data[]
        label: "normal_sequence: WT sequencing data and readgroup information"
        doc: |
          normal_sequence represents the sequencing data for the WT sample as either FASTQs or BAMs with
          accompanying readgroup information. Note that in the @RG field ID and SM are required.
    normal_name:
        type: string?
        default: 'normal'
    trimming:
        type:
            - ../types/trimming_options.yml#trimming_options
            - "null"
    bqsr_known_sites:
        type: File[]
        secondaryFiles: [.tbi]
        doc: "One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis."
    bqsr_intervals:
        type: string[]
    target_intervals:
        type: File
    per_base_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    per_target_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    qc_intervals:
        type: File
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
    strelka_cpu_reserved:
        type: int?
        default: 8
    scatter_count:
        type: int
        doc: "scatters each supported variant detector (varscan, pindel, mutect) into this many parallel jobs"
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
    filter_somatic_llr_threshold:
        type: float
        default: 5
        doc: "Sets the stringency (log-likelihood ratio) used to filter out non-somatic variants.  Typical values are 10=high stringency, 5=normal, 3=low stringency. Low stringency may be desirable when read depths are low (as in WGS) or when tumor samples are impure."
    filter_somatic_llr_tumor_purity:
        type: float
        default: 1
        doc: "Sets the purity of the tumor used in the somatic llr filter, used to remove non-somatic variants. Probably only needs to be adjusted for low-purity (< 50%).  Range is 0 to 1"
    filter_somatic_llr_normal_contamination_rate:
        type: float
        default: 0
        doc: "Sets the fraction of tumor present in the normal sample (range 0 to 1), used in the somatic llr filter. Useful for heavily contaminated adjacent normals. Range is 0 to 1"
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
    vep_custom_annotations:
        type: ../types/vep_custom_annotation.yml#vep_custom_annotation[]
        doc: "custom type, check types directory for input format"
    manta_call_regions:
        type: File?
        secondaryFiles: [.tbi]
    manta_non_wgs:
        type: boolean?
        default: false
    manta_output_contigs:
        type: boolean?
    somalier_vcf:
        type: File
    tumor_sample_name:
        type: string
    normal_sample_name:
        type: string
    validated_variants:
        type: File?
        secondaryFiles: [.tbi]
        doc: "An optional VCF with variants that will be flagged as 'VALIDATED' if found in this pipeline's main output VCF"
    cnvkit_target_average_size:
        type: int?
        doc: "approximate size of split target bins for CNVkit; if not set a suitable window size will be set by CNVkit automatically"
outputs:
##tumor alignment and QC
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
    tumor_per_target_coverage_metrics:
        type: File[]
        outputSource: tumor_alignment_and_qc/per_target_coverage_metrics
    tumor_per_target_hs_metrics:
        type: File[]
        outputSource: tumor_alignment_and_qc/per_target_hs_metrics
    tumor_per_base_coverage_metrics:
        type: File[]
        outputSource: tumor_alignment_and_qc/per_base_coverage_metrics
    tumor_per_base_hs_metrics:
        type: File[]
        outputSource: tumor_alignment_and_qc/per_base_hs_metrics
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
    tumor_insert_size_histogram:
        type: File
        outputSource: tumor_alignment_and_qc/insert_size_histogram
    tumor_gc_bias_metrics:
        type: File
        outputSource: tumor_alignment_and_qc/gc_bias_metrics
    tumor_gc_bias_metrics_chart:
        type: File
        outputSource: tumor_alignment_and_qc/gc_bias_metrics_chart
    tumor_gc_bias_metrics_summary:
        type: File
        outputSource: tumor_alignment_and_qc/gc_bias_metrics_summary
    tumor_wgs_metrics:
        type: File
        outputSource: tumor_alignment_and_qc/wgs_metrics
##normal alignment and qc
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
    normal_per_target_coverage_metrics:
        type: File[]
        outputSource: normal_alignment_and_qc/per_target_coverage_metrics
    normal_per_target_hs_metrics:
        type: File[]
        outputSource: normal_alignment_and_qc/per_target_hs_metrics
    normal_per_base_coverage_metrics:
        type: File[]
        outputSource: normal_alignment_and_qc/per_base_coverage_metrics
    normal_per_base_hs_metrics:
        type: File[]
        outputSource: normal_alignment_and_qc/per_base_hs_metrics
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
    normal_insert_size_histogram:
        type: File
        outputSource: normal_alignment_and_qc/insert_size_histogram
    normal_gc_bias_metrics:
        type: File
        outputSource: normal_alignment_and_qc/gc_bias_metrics
    normal_gc_bias_metrics_chart:
        type: File
        outputSource: normal_alignment_and_qc/gc_bias_metrics_chart
    normal_gc_bias_metrics_summary:
        type: File
        outputSource: normal_alignment_and_qc/gc_bias_metrics_summary
    normal_wgs_metrics:
        type: File
        outputSource: normal_alignment_and_qc/wgs_metrics
##variant calling
    mutect_unfiltered_vcf:
        type: File
        outputSource: detect_variants/mutect_unfiltered_vcf
        secondaryFiles: [.tbi]
    mutect_filtered_vcf:
        type: File
        outputSource: detect_variants/mutect_filtered_vcf
        secondaryFiles: [.tbi]
    strelka_unfiltered_vcf:
        type: File
        outputSource: detect_variants/strelka_unfiltered_vcf
        secondaryFiles: [.tbi]
    strelka_filtered_vcf:
        type: File
        outputSource: detect_variants/strelka_filtered_vcf
        secondaryFiles: [.tbi]
    varscan_unfiltered_vcf:
        type: File
        outputSource: detect_variants/varscan_unfiltered_vcf
        secondaryFiles: [.tbi]
    varscan_filtered_vcf:
        type: File
        outputSource: detect_variants/varscan_filtered_vcf
        secondaryFiles: [.tbi]
    docm_filtered_vcf:
        type: File
        outputSource: detect_variants/docm_filtered_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: detect_variants/final_vcf
        secondaryFiles: [.tbi]
    final_filtered_vcf:
        type: File
        outputSource: detect_variants/final_filtered_vcf
        secondaryFiles: [.tbi]
    final_tsv:
        type: File
        outputSource: detect_variants/final_tsv
    vep_summary:
        type: File
        outputSource: detect_variants/vep_summary
    tumor_snv_bam_readcount_tsv:
        type: File
        outputSource: detect_variants/tumor_snv_bam_readcount_tsv
    tumor_indel_bam_readcount_tsv:
        type: File
        outputSource: detect_variants/tumor_indel_bam_readcount_tsv
    normal_snv_bam_readcount_tsv:
        type: File
        outputSource: detect_variants/normal_snv_bam_readcount_tsv
    normal_indel_bam_readcount_tsv:
        type: File
        outputSource: detect_variants/normal_indel_bam_readcount_tsv
    diploid_variants:
        type: File?
        outputSource: manta/diploid_variants
        secondaryFiles: [.tbi]
    somatic_variants:
        type: File?
        outputSource: manta/somatic_variants
        secondaryFiles: [.tbi]
    all_candidates:
        type: File
        outputSource: manta/all_candidates
        secondaryFiles: [.tbi]
    small_candidates:
        type: File
        outputSource: manta/small_candidates
        secondaryFiles: [.tbi]
    tumor_only_variants:
        type: File?
        outputSource: manta/tumor_only_variants
        secondaryFiles: [.tbi]
##sample concordance check
    somalier_concordance_metrics:
        type: File
        outputSource: concordance/somalier_pairs
    somalier_concordance_statistics:
        type: File
        outputSource: concordance/somalier_samples
steps:
    tumor_alignment_and_qc:
        run: alignment_wgs.cwl
        in:
            reference: reference
            sequence: tumor_sequence
            trimming: trimming
            omni_vcf: omni_vcf
            intervals: qc_intervals
            picard_metric_accumulation_level: picard_metric_accumulation_level
            bqsr_known_sites: bqsr_known_sites
            bqsr_intervals: bqsr_intervals
            minimum_mapping_quality: qc_minimum_mapping_quality
            minimum_base_quality: qc_minimum_base_quality
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
            sample_name: tumor_name
        out:
            [alignment_summary_metrics, bam, flagstats, gc_bias_metrics_chart, gc_bias_metrics_summary, gc_bias_metrics, insert_size_histogram, insert_size_metrics, mark_duplicates_metrics, per_base_coverage_metrics, per_base_hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, summary_hs_metrics, verify_bam_id_depth, verify_bam_id_metrics, wgs_metrics]
    normal_alignment_and_qc:
        run: alignment_wgs.cwl
        in:
            reference: reference
            sequence: normal_sequence
            trimming: trimming
            omni_vcf: omni_vcf
            intervals: qc_intervals
            picard_metric_accumulation_level: picard_metric_accumulation_level
            bqsr_known_sites: bqsr_known_sites
            bqsr_intervals: bqsr_intervals
            minimum_mapping_quality: qc_minimum_mapping_quality
            minimum_base_quality: qc_minimum_base_quality
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
            sample_name: normal_name
        out:
            [alignment_summary_metrics, bam, flagstats, gc_bias_metrics_chart, gc_bias_metrics_summary, gc_bias_metrics, insert_size_histogram, insert_size_metrics, mark_duplicates_metrics, per_base_coverage_metrics, per_base_hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, summary_hs_metrics, verify_bam_id_depth, verify_bam_id_metrics, wgs_metrics]
    concordance:
        run: ../tools/concordance.cwl
        in:
            reference: reference
            bam_1: tumor_alignment_and_qc/bam
            bam_2: normal_alignment_and_qc/bam
            vcf: somalier_vcf
        out:
            [somalier_pairs, somalier_samples]
    detect_variants:
        run: detect_variants_wgs.cwl
        in:
            reference: reference
            tumor_bam: tumor_alignment_and_qc/bam
            normal_bam: normal_alignment_and_qc/bam
            roi_intervals: target_intervals
            strelka_exome_mode:
                default: false
            strelka_cpu_reserved: strelka_cpu_reserved
            scatter_count: scatter_count
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
            filter_somatic_llr_threshold: filter_somatic_llr_threshold
            filter_somatic_llr_tumor_purity: filter_somatic_llr_tumor_purity
            filter_somatic_llr_normal_contamination_rate: filter_somatic_llr_normal_contamination_rate
            vep_cache_dir: vep_cache_dir
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            vep_pick: vep_pick
            cle_vcf_filter: cle_vcf_filter
            variants_to_table_fields: variants_to_table_fields
            variants_to_table_genotype_fields: variants_to_table_genotype_fields
            vep_to_table_fields: vep_to_table_fields
            tumor_sample_name: tumor_sample_name
            normal_sample_name: normal_sample_name
            vep_custom_annotations: vep_custom_annotations
            validated_variants: validated_variants
        out:
            [mutect_unfiltered_vcf, mutect_filtered_vcf, strelka_unfiltered_vcf, strelka_filtered_vcf, varscan_unfiltered_vcf, varscan_filtered_vcf, docm_filtered_vcf, final_vcf, final_filtered_vcf, final_tsv, vep_summary, tumor_snv_bam_readcount_tsv, tumor_indel_bam_readcount_tsv, normal_snv_bam_readcount_tsv, normal_indel_bam_readcount_tsv]
    cnvkit:
        run: ../tools/cnvkit_batch.cwl
        in:
            tumor_bam: tumor_alignment_and_qc/bam
            method:
                default: 'wgs'
            reference:
                source: [normal_alignment_and_qc/bam, reference]
                valueFrom: |
                    ${
                      var normal = self[0];
                      var fasta = self[1];
                      return {'normal_bam': normal, 'fasta_file': fasta};
                    }
            target_average_size: cnvkit_target_average_size
        out:
            [intervals_antitarget, intervals_target, normal_antitarget_coverage, normal_target_coverage, reference_coverage, cn_diagram, cn_scatter_plot, tumor_antitarget_coverage, tumor_target_coverage, tumor_bin_level_ratios, tumor_segmented_ratios]
    manta: 
        run: ../tools/manta_somatic.cwl
        in:
            normal_bam: normal_alignment_and_qc/bam
            tumor_bam: tumor_alignment_and_qc/bam
            reference: reference
            non_wgs: manta_non_wgs
            output_contigs: manta_output_contigs
        out:
            [diploid_variants, somatic_variants, all_candidates, small_candidates, tumor_only_variants]
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

