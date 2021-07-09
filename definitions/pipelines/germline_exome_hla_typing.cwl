#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "exome alignment and germline variant detection, with optitype for HLA typing"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
          - $import: ../types/sequence_data.yml
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
        label: "sequence: sequencing data and readgroup information"
        doc: |
          sequence represents the sequencing data as either FASTQs or BAMs with accompanying
          readgroup information. Note that in the @RG field ID and SM are required.
    bqsr_known_sites:
        type: File[]
        secondaryFiles: [.tbi]
        doc: "One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis."
    bqsr_intervals:
        type: string[]?
    bait_intervals:
        type: File
        label: "bait_intervals: interval_list file of baits used in the sequencing experiment"
        doc: |
          bait_intervals is an interval_list corresponding to the baits used in sequencing reagent.
          These are essentially coordinates for regions you were able to design probes for in the reagent.
          Typically the reagent provider has this information available in bed format and it can be
          converted to an interval_list with Picard BedToIntervalList. Astrazeneca also maintains a repo
          of baits for common sequencing reagents available at https://github.com/AstraZeneca-NGS/reference_data
    target_intervals:
        type: File
        label: "target_intervals: interval_list file of targets used in the sequencing experiment"
        doc: |
          target_intervals is an interval_list corresponding to the targets for the capture reagent.
          Bed files with this information can be converted to interval_lists with Picard BedToIntervalList.
          In general for a WES exome reagent bait_intervals and target_intervals are the same.
    target_interval_padding:
        type: int
        label: "target_interval_padding: number of bp flanking each target region in which to allow variant calls"
        doc: |
            The effective coverage of capture products generally extends out beyond the actual regions
            targeted. This parameter allows variants to be called in these wingspan regions, extending
            this many base pairs from each side of the target regions.
        default: 100
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
    emit_reference_confidence:
        type:
            type: enum
            symbols: ['NONE', 'BP_RESOLUTION', 'GVCF']
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
    vep_custom_annotations:
        type: ../types/vep_custom_annotation.yml#vep_custom_annotation[]
        doc: "custom type, check types directory for input format"
    qc_minimum_mapping_quality:
        type: int?
    qc_minimum_base_quality:
        type: int?
    optitype_name:
        type: string?
outputs:
    cram:
        type: File
        outputSource: germline_exome/cram
    mark_duplicates_metrics:
        type: File
        outputSource: germline_exome/mark_duplicates_metrics
    insert_size_metrics:
        type: File
        outputSource: germline_exome/insert_size_metrics
    insert_size_histogram:
        type: File
        outputSource: germline_exome/insert_size_histogram
    alignment_summary_metrics:
        type: File
        outputSource: germline_exome/alignment_summary_metrics
    hs_metrics:
        type: File
        outputSource: germline_exome/hs_metrics
    per_target_coverage_metrics:
        type: File[]
        outputSource: germline_exome/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File[]
        outputSource: germline_exome/per_target_hs_metrics
    per_base_coverage_metrics:
        type: File[]
        outputSource: germline_exome/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File[]
        outputSource: germline_exome/per_base_hs_metrics
    summary_hs_metrics:
        type: File[]
        outputSource: germline_exome/summary_hs_metrics
    flagstats:
        type: File
        outputSource: germline_exome/flagstats
    verify_bam_id_metrics:
        type: File
        outputSource: germline_exome/verify_bam_id_metrics
    verify_bam_id_depth:
        type: File
        outputSource: germline_exome/verify_bam_id_depth
    raw_vcf:
        type: File
        outputSource: germline_exome/raw_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: germline_exome/final_vcf
        secondaryFiles: [.tbi]
    filtered_vcf:
        type: File
        outputSource: germline_exome/filtered_vcf
        secondaryFiles: [.tbi]
    vep_summary:
        type: File
        outputSource: germline_exome/vep_summary
    optitype_tsv:
        type: File
        outputSource: optitype/optitype_tsv
    optitype_plot:
        type: File
        outputSource: optitype/optitype_plot
steps:
    germline_exome:
        run: germline_exome.cwl
        in:
            reference: reference
            sequence: sequence
            bqsr_known_sites: bqsr_known_sites
            bqsr_intervals: bqsr_intervals
            bait_intervals: bait_intervals
            target_intervals: target_intervals
            target_interval_padding: target_interval_padding
            per_base_intervals: per_base_intervals
            per_target_intervals: per_target_intervals
            summary_intervals: summary_intervals
            omni_vcf: omni_vcf
            picard_metric_accumulation_level: picard_metric_accumulation_level
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            ploidy: ploidy
            vep_cache_dir: vep_cache_dir
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            vep_custom_annotations: vep_custom_annotations
            qc_minimum_mapping_quality: qc_minimum_mapping_quality
            qc_minimum_base_quality: qc_minimum_base_quality
        out:
            [cram, mark_duplicates_metrics, insert_size_metrics, insert_size_histogram, alignment_summary_metrics, hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, per_base_coverage_metrics, per_base_hs_metrics, summary_hs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth, raw_vcf, final_vcf, filtered_vcf, vep_summary]
    optitype:
        run: ../tools/optitype_dna.cwl
        in:
            optitype_name: optitype_name
            cram: germline_exome/cram
            reference: reference
        out:
            [optitype_tsv, optitype_plot]
