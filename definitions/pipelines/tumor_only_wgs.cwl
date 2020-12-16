#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "wgs alignment and tumor-only variant detection"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
          - $import: ../types/sequence_data.yml
          - $import: ../types/trimming_options.yml
          - $import: ../types/vep_custom_annotation.yml
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
    sequence:
        type: ../types/sequence_data.yml#sequence_data[]
    trimming:
        type:
            - ../types/trimming_options.yml#trimming_options
            - "null"
    bqsr_known_sites:
        type: File[]
        secondaryFiles: [.tbi]
        doc: "One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis."
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
    bqsr_intervals:
        type: string[]
        label: "bqsr_intervals: Array of strings specifying regions for base quality score recalibration"
        doc: |
          bqsr_intervals provides an array of genomic intervals for which to apply
          GATK base quality score recalibrations. Typically intervals are given
          for the entire chromosome (i.e. chr1, chr2, etc.), these names should match
          the format in the reference file.
    target_intervals:
        type: File
        label: "target_intervals: interval_list file of targets used in the sequencing experiment"
        doc: |
          target_intervals is an interval_list corresponding to the targets for the reagent. In the case of WGS
          this generally should contain intervals that span the entire genome. 
    per_base_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
        label: "per_base_intervals: additional intervals over which to summarize coverage/QC at a per-base resolution"
        doc: "per_base_intervals is a list of regions (in interval_list format) over which to summarize coverage/QC at a per-base resolution."
    per_target_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
        label: "per_target_intervals: additional intervals over which to summarize coverage/QC at a per-target resolution"
        doc: "per_target_intervals list of regions (in interval_list format) over which to summarize coverage/QC at a per-target resolution."
    summary_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    picard_metric_accumulation_level:
        type: string
    roi_intervals:
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
    synonyms_file:
        type: File?
    vep_pick:
        type:
            - "null"
            - type: enum
              symbols: ["pick", "flag_pick", "pick_allele", "per_gene", "pick_allele_gene", "flag_pick_allele", "flag_pick_allele_gene"]
    sample_name:
        type: string
    docm_vcf:
        type: File
        doc: "Common mutations in cancer that will be genotyped and passed through into the merged VCF if they have even low-level evidence of a mutation (by default, marked with filter DOCM_ONLY)"
    vep_custom_annotations:
        type: ../types/vep_custom_annotation.yml#vep_custom_annotation[]
        doc: "custom type, check types directory for input format"
    readcount_minimum_mapping_quality:
        type: int?
    readcount_minimum_base_quality:
        type: int?
    variants_to_table_fields:
        type: string[]
        default: [CHROM,POS,ID,REF,ALT,set,AC,AF]
    variants_to_table_genotype_fields:
        type: string[]
        default: [GT,AD]
    vep_to_table_fields:
        type: string[]
        default: [HGVSc,HGVSp]
    varscan_min_coverage:
        type: int?
        default: 8
    varscan_min_var_freq:
        type: float?
        default: 0.05
    varscan_min_reads:
        type: int?
        default: 2
    maximum_population_allele_frequency:
        type: float?
        default: 0.001
    qc_minimum_mapping_quality:
        type: int?
        default: 0
    qc_minimum_base_quality:
        type: int?
        default: 0
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
    varscan_vcf:
        type: File
        outputSource: detect_variants/varscan_vcf
        secondaryFiles: [.tbi]
    docm_gatk_vcf:
        type: File
        outputSource: detect_variants/docm_gatk_vcf
    annotated_vcf:
        type: File
        outputSource: detect_variants/annotated_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: detect_variants/final_vcf
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
steps:
    alignment_and_qc:
        run: alignment_wgs.cwl
        in:
            reference: reference
            sequence: sequence
            trimming: trimming
            omni_vcf: omni_vcf
            intervals: target_intervals
            picard_metric_accumulation_level: picard_metric_accumulation_level
            bqsr_known_sites: bqsr_known_sites
            bqsr_intervals: bqsr_intervals
            minimum_mapping_quality: qc_minimum_mapping_quality
            minimum_base_quality: qc_minimum_base_quality
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
            sample_name: sample_name
        out:
            [bam, mark_duplicates_metrics, insert_size_metrics, insert_size_histogram, alignment_summary_metrics, gc_bias_metrics, gc_bias_metrics_chart, gc_bias_metrics_summary, wgs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth, per_base_coverage_metrics, per_base_hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, summary_hs_metrics]
    detect_variants:
        run: tumor_only_detect_variants.cwl
        in:
            reference: reference
            bam: alignment_and_qc/bam
            roi_intervals: roi_intervals
            varscan_min_coverage: varscan_min_coverage
            varscan_min_var_freq: varscan_min_var_freq
            varscan_min_reads: varscan_min_reads
            maximum_population_allele_frequency: maximum_population_allele_frequency
            vep_cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            vep_pick: vep_pick
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            variants_to_table_fields: variants_to_table_fields
            variants_to_table_genotype_fields: variants_to_table_genotype_fields
            vep_to_table_fields: vep_to_table_fields
            sample_name: sample_name
            docm_vcf: docm_vcf
            vep_custom_annotations: vep_custom_annotations
            readcount_minimum_mapping_quality: readcount_minimum_mapping_quality
            readcount_minimum_base_quality: readcount_minimum_base_quality
        out:
            [varscan_vcf, docm_gatk_vcf, annotated_vcf, final_vcf, final_tsv, vep_summary, tumor_snv_bam_readcount_tsv, tumor_indel_bam_readcount_tsv]
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
