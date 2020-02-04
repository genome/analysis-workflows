#/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Immunotherapy Workflow"
requirements:
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
          - $import: ../types/sequence_data.yml
          - $import: ../types/vep_custom_annotation.yml
    - class: SubworkflowFeatureRequirement
inputs:
    #rnaseq inputs
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
        type:
          - "null"
          - type: enum
            symbols: ["first", "second", "unstranded"]
    refFlat:
        type: File
    ribosomal_intervals:
        type: File?

    #somatic inputs
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict, .amb, .ann, .bwt, .pac, .sa]
        label: "reference: Reference fasta file for a desired assembly"
        doc: |
          reference contains the nucleotide sequence for a given assembly (hg37, hg38, etc.)
          in fasta format for the entire genome. This is what reads will be aligned to.
          Appropriate files can be found on ensembl at https://ensembl.org/info/data/ftp/index.html
          When providing the reference secondary files corresponding to reference indices must be
          located in the same directory as the reference itself. These files can be created with
          samtools index, bwa index, and picard CreateSequenceDictionary.

    tumor_sequence:
        type: ../types/sequence_data.yml#sequence_data[]
        label: "tumor_sequence: file specifying the location of MT sequencing data"
        doc: |
          tumor_sequence is a data structure described in sequence_data.yml used to pass information regarding
          sequencing data for single sample (i.e. fastq files). If more than one fastq file exist
          for a sample, as in the case for multiple instrument data, the sequence tag is simply
          repeated with the additional data (see example input file). Note that in the @RG field
          ID and SM are required.
    tumor_name:
        type: string?
        default: 'tumor'
        label: "tumor_name: String specifying the name of the MT sample"
        doc: |
          tumor_name provides a string for what the MT sample will be referred to in the various
          outputs, for exmaple the VCF files.
    normal_sequence:
        type: ../types/sequence_data.yml#sequence_data[]
        label: "normal_sequence: file specifying the location of WT sequencing data"
        doc: |
          normal_sequence is a data structure described in sequence_data.yml used to pass information regarding
          sequencing data for single sample (i.e. fastq files). If more than one fastq file exist
          for a sample, as in the case for multiple instrument data, the sequence tag is simply
          repeated with the additional data (see example input file). Note that in the @RG field
          ID and SM are required.
    normal_name:
        type: string?
        default: 'normal'
        label: "normal_name: String specifying the name of the WT sample"
        doc: |
          normal_name provides a string for what the WT sample will be referred to in the various
          outputs, for exmaple the VCF files.
    mills:
        type: File
        secondaryFiles: [.tbi]
        label: "mills: File specifying common polymorphic indels from mills et al."
        doc: |
          mills provides known polymorphic indels recommended by GATK for a variety of
          tools including the BaseRecalibrator. This file is part of the GATK resource
          bundle available at http://www.broadinstitute.org/gatk/guide/article?id=1213
          Essentially it is a list of known indels originally discovered by mill et al.
          https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1557762/
          File should be in vcf format, and tabix indexed.
    known_indels:
        type: File
        secondaryFiles: [.tbi]
        label: "known_indels: File specifying common polymorphic indels from 1000G"
        doc: |
          known_indels provides known indels reecommended by GATK for a variety of tools
          including the BaseRecalibrator. This file is part of the GATK resource bundle
          available at http://www.broadinstitute.org/gatk/guide/article?id=1213
          Essintially it is a list of known indels from 1000 Genomes Phase I indel calls.
          File should be in vcf format, and tabix indexed.
    dbsnp_vcf:
        type: File
        secondaryFiles: [.tbi]
        label: "dbsnp_vcf: File specifying common polymorphic indels from dbSNP"
        doc: |
          dbsnp_vcf provides known indels reecommended by GATK for a variety of tools
          including the BaseRecalibrator. This file is part of the GATK resource bundle
          available at http://www.broadinstitute.org/gatk/guide/article?id=1213
          Essintially it is a list of known indels from dbSNP. File should be in vcf format,
          and tabix indexed.
    bqsr_intervals:
        type: string[]
        label: "bqsr_intervals: Array of strings specifying regions for base quality score recalibration"
        doc: |
          bqsr_intervals provides an array of genomic intervals for which to apply
          GATK base quality score recalibrations. Typically intervals are given
          for the entire chromosome (i.e. chr1, chr2, etc.), these names should match
          the format in the reference file.
    bait_intervals:
        type: File
        label: "bait_intervals: interval_list file of baits used in the sequencing experiment"
        doc: |
          bait_intervals is an interval_list corresponding to the baits used in sequencing reagent.
          These are essentially coordinates for regions you were able to design probes for in the reagent.
          Typically the reagent provider has this information available in bed format and it can be
          converted to an interval_list with Picards BedToIntervalList. AstraZeneca also maintains a repo
          of baits for common sequencing reagents available at https://github.com/AstraZeneca-NGS/reference_data
    target_intervals:
        type: File
        label: "target_intervals: interval_list file of targets used in the sequencing experiment"
        doc: |
          target_intervals is an interval_list corresponding to the targets for the sequencing reagent.
          These are essentially coordinates for regions you wanted to design probes for in the reagent.
          Bed files with this information can be converted to interval_lists with Picards BedToIntervalList.
          In general for a WES exome reagent bait_intervals and target_intervals are the same.
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
        default: true
    manta_output_contigs:
        type: boolean?
    somalier_vcf:
        type: File
    known_variants:
        type: File?
        secondaryFiles: [.tbi]
        doc: "Previously discovered variants to be flagged in this pipelines's output vcf"

    #germline inputs
    emit_reference_confidence:
        type:
            type: enum
            symbols: ['NONE', 'BP_RESOLUTION', 'GVCF']
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

    #phase_vcf inputs
    reference_dict:
        type: File

    clinical_calls:
        type: string[]?
        doc: "Clinical HLA typing results; element format: HLA-X*01:02[/HLA-X...]"

    #pvacseq inputs
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

    tumor_sample_name:
        type: string
    normal_sample_name:
        type: string

outputs:
    final_bam:
        type: File
        outputSource: rnaseq/final_bam
        secondaryFiles: [.bai]
    stringtie_transcript_gtf:
        type: File
        outputSource: rnaseq/stringtie_transcript_gtf
    stringtie_gene_expression_tsv:
        type: File
        outputSource: rnaseq/stringtie_gene_expression_tsv
    transcript_abundance_tsv:
        type: File
        outputSource: rnaseq/transcript_abundance_tsv
    transcript_abundance_h5:
        type: File
        outputSource: rnaseq/transcript_abundance_h5
    gene_abundance:
        type: File
        outputSource: rnaseq/gene_abundance
    metrics:
        type: File
        outputSource: rnaseq/metrics
    chart:
        type: File?
        outputSource: rnaseq/chart

    tumor_cram:
        type: File
        outputSource: somatic/tumor_cram
    tumor_mark_duplicates_metrics:
        type: File
        outputSource: somatic/tumor_mark_duplicates_metrics
    tumor_insert_size_metrics:
        type: File
        outputSource: somatic/tumor_insert_size_metrics
    tumor_alignment_summary_metrics:
        type: File
        outputSource: somatic/tumor_alignment_summary_metrics
    tumor_hs_metrics:
        type: File
        outputSource: somatic/tumor_hs_metrics
    tumor_per_target_coverage_metrics:
        type: File[]
        outputSource: somatic/tumor_per_target_coverage_metrics
    tumor_per_target_hs_metrics:
        type: File[]
        outputSource: somatic/tumor_per_target_hs_metrics
    tumor_per_base_coverage_metrics:
        type: File[]
        outputSource: somatic/tumor_per_base_coverage_metrics
    tumor_per_base_hs_metrics:
        type: File[]
        outputSource: somatic/tumor_per_base_hs_metrics
    tumor_summary_hs_metrics:
        type: File[]
        outputSource: somatic/tumor_summary_hs_metrics
    tumor_flagstats:
        type: File
        outputSource: somatic/tumor_flagstats
    tumor_verify_bam_id_metrics:
        type: File
        outputSource: somatic/tumor_verify_bam_id_metrics
    tumor_verify_bam_id_depth:
        type: File
        outputSource: somatic/tumor_verify_bam_id_depth
    normal_cram:
        type: File
        outputSource: somatic/normal_cram
    normal_mark_duplicates_metrics:
        type: File
        outputSource: somatic/normal_mark_duplicates_metrics
    normal_insert_size_metrics:
        type: File
        outputSource: somatic/normal_insert_size_metrics
    normal_alignment_summary_metrics:
        type: File
        outputSource: somatic/normal_alignment_summary_metrics
    normal_hs_metrics:
        type: File
        outputSource: somatic/normal_hs_metrics
    normal_per_target_coverage_metrics:
        type: File[]
        outputSource: somatic/normal_per_target_coverage_metrics
    normal_per_target_hs_metrics:
        type: File[]
        outputSource: somatic/normal_per_target_hs_metrics
    normal_per_base_coverage_metrics:
        type: File[]
        outputSource: somatic/normal_per_base_coverage_metrics
    normal_per_base_hs_metrics:
        type: File[]
        outputSource: somatic/normal_per_base_hs_metrics
    normal_summary_hs_metrics:
        type: File[]
        outputSource: somatic/normal_summary_hs_metrics
    normal_flagstats:
        type: File
        outputSource: somatic/normal_flagstats
    normal_verify_bam_id_metrics:
        type: File
        outputSource: somatic/normal_verify_bam_id_metrics
    normal_verify_bam_id_depth:
        type: File
        outputSource: somatic/normal_verify_bam_id_depth
    mutect_unfiltered_vcf:
        type: File
        outputSource: somatic/mutect_unfiltered_vcf
        secondaryFiles: [.tbi]
    mutect_filtered_vcf:
        type: File
        outputSource: somatic/mutect_filtered_vcf
        secondaryFiles: [.tbi]
    strelka_unfiltered_vcf:
        type: File
        outputSource: somatic/strelka_unfiltered_vcf
        secondaryFiles: [.tbi]
    strelka_filtered_vcf:
        type: File
        outputSource: somatic/strelka_filtered_vcf
        secondaryFiles: [.tbi]
    varscan_unfiltered_vcf:
        type: File
        outputSource: somatic/varscan_unfiltered_vcf
        secondaryFiles: [.tbi]
    varscan_filtered_vcf:
        type: File
        outputSource: somatic/varscan_filtered_vcf
        secondaryFiles: [.tbi]
    pindel_unfiltered_vcf:
        type: File
        outputSource: somatic/pindel_unfiltered_vcf
        secondaryFiles: [.tbi]
    pindel_filtered_vcf:
        type: File
        outputSource: somatic/pindel_filtered_vcf
        secondaryFiles: [.tbi]
    docm_filtered_vcf:
        type: File
        outputSource: somatic/docm_filtered_vcf
        secondaryFiles: [.tbi]
    somatic_final_vcf:
        type: File
        outputSource: somatic/final_vcf
        secondaryFiles: [.tbi]
    final_filtered_vcf:
        type: File
        outputSource: somatic/final_filtered_vcf
        secondaryFiles: [.tbi]
    final_tsv:
        type: File
        outputSource: somatic/final_tsv
    somatic_vep_summary:
        type: File
        outputSource: somatic/vep_summary
    tumor_snv_bam_readcount_tsv:
        type: File
        outputSource: somatic/tumor_snv_bam_readcount_tsv
    tumor_indel_bam_readcount_tsv:
        type: File
        outputSource: somatic/tumor_indel_bam_readcount_tsv
    normal_snv_bam_readcount_tsv:
        type: File
        outputSource: somatic/normal_snv_bam_readcount_tsv
    normal_indel_bam_readcount_tsv:
        type: File
        outputSource: somatic/normal_indel_bam_readcount_tsv
    intervals_antitarget:
        type: File?
        outputSource: somatic/intervals_antitarget
    intervals_target:
        type: File?
        outputSource: somatic/intervals_target
    normal_antitarget_coverage:
        type: File
        outputSource: somatic/normal_antitarget_coverage
    normal_target_coverage:
        type: File
        outputSource: somatic/normal_target_coverage
    reference_coverage:
        type: File?
        outputSource: somatic/reference_coverage
    cn_diagram:
        type: File?
        outputSource: somatic/cn_diagram
    cn_scatter_plot:
        type: File?
        outputSource: somatic/cn_scatter_plot
    tumor_antitarget_coverage:
        type: File
        outputSource: somatic/tumor_antitarget_coverage
    tumor_target_coverage:
        type: File
        outputSource: somatic/tumor_target_coverage
    tumor_bin_level_ratios:
        type: File
        outputSource: somatic/tumor_bin_level_ratios
    tumor_segmented_ratios:
        type: File
        outputSource: somatic/tumor_segmented_ratios
    diploid_variants:
        type: File?
        outputSource: somatic/diploid_variants
        secondaryFiles: [.tbi]
    somatic_variants:
        type: File?
        outputSource: somatic/somatic_variants
        secondaryFiles: [.tbi]
    all_candidates:
        type: File
        outputSource: somatic/all_candidates
        secondaryFiles: [.tbi]
    small_candidates:
        type: File
        outputSource: somatic/small_candidates
        secondaryFiles: [.tbi]
    tumor_only_variants:
        type: File?
        outputSource: somatic/tumor_only_variants
        secondaryFiles: [.tbi]
    somalier_concordance_metrics:
        type: File
        outputSource: somatic/somalier_concordance_metrics
    somalier_concordance_statistics:
        type: File
        outputSource: somatic/somalier_concordance_statistics

    cram:
        type: File
        outputSource: germline/cram
    mark_duplicates_metrics:
        type: File
        outputSource: germline/mark_duplicates_metrics
    insert_size_metrics:
        type: File
        outputSource: germline/insert_size_metrics
    insert_size_histogram:
        type: File
        outputSource: germline/insert_size_histogram
    alignment_summary_metrics:
        type: File
        outputSource: germline/alignment_summary_metrics
    hs_metrics:
        type: File
        outputSource: germline/hs_metrics
    per_target_coverage_metrics:
        type: File[]
        outputSource: germline/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File[]
        outputSource: germline/per_target_hs_metrics
    per_base_coverage_metrics:
        type: File[]
        outputSource: germline/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File[]
        outputSource: germline/per_base_hs_metrics
    summary_hs_metrics:
        type: File[]
        outputSource: germline/summary_hs_metrics
    flagstats:
        type: File
        outputSource: germline/flagstats
    verify_bam_id_metrics:
        type: File
        outputSource: germline/verify_bam_id_metrics
    verify_bam_id_depth:
        type: File
        outputSource: germline/verify_bam_id_depth
    gvcf:
        type: File[]
        outputSource: germline/gvcf
    germline_final_vcf:
        type: File
        outputSource: germline/final_vcf
        secondaryFiles: [.tbi]
    germline_filtered_vcf:
        type: File
        outputSource: germline/filtered_vcf
        secondaryFiles: [.tbi]
    germline_vep_summary:
        type: File
        outputSource: germline/vep_summary
    optitype_tsv:
        type: File
        outputSource: germline/optitype_tsv
    optitype_plot:
        type: File
        outputSource: germline/optitype_plot

    phased_vcf:
        type: File
        outputSource: phase_vcf/phased_vcf
        secondaryFiles: [.tbi]

    allele_string:
        type: string[]
        outputSource: extract_alleles/allele_string
    consensus_alleles:
        type: string[]
        outputSource: hla_consensus/consensus_alleles
    hla_call_files:
        type: Directory
        outputSource: hla_consensus/hla_call_files

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
            species: vep_ensembl_species
            assembly: vep_ensembl_assembly
        out:
            [final_bam, stringtie_transcript_gtf, stringtie_gene_expression_tsv, transcript_abundance_tsv, transcript_abundance_h5, gene_abundance, metrics, chart, fusion_evidence]
    somatic:
        run: somatic_exome.cwl
        in:
            reference: reference
            tumor_sequence: tumor_sequence
            tumor_name: tumor_name
            normal_sequence: normal_sequence
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
            vep_custom_annotations: vep_custom_annotations
            manta_call_regions: manta_call_regions
            manta_non_wgs: manta_non_wgs
            manta_output_contigs: manta_output_contigs
            somalier_vcf: somalier_vcf
            tumor_sample_name: tumor_sample_name
            normal_sample_name: normal_sample_name
            known_variants: known_variants
        out:
            [tumor_cram,tumor_mark_duplicates_metrics,tumor_insert_size_metrics,tumor_alignment_summary_metrics,tumor_hs_metrics,tumor_per_target_coverage_metrics,tumor_per_target_hs_metrics,tumor_per_base_coverage_metrics,tumor_per_base_hs_metrics,tumor_summary_hs_metrics,tumor_flagstats,tumor_verify_bam_id_metrics,tumor_verify_bam_id_depth,normal_cram,normal_mark_duplicates_metrics,normal_insert_size_metrics,normal_alignment_summary_metrics,normal_hs_metrics,normal_per_target_coverage_metrics,normal_per_target_hs_metrics,normal_per_base_coverage_metrics,normal_per_base_hs_metrics,normal_summary_hs_metrics,normal_flagstats,normal_verify_bam_id_metrics,normal_verify_bam_id_depth,mutect_unfiltered_vcf,mutect_filtered_vcf,strelka_unfiltered_vcf,strelka_filtered_vcf,varscan_unfiltered_vcf,varscan_filtered_vcf,pindel_unfiltered_vcf,pindel_filtered_vcf,docm_filtered_vcf,final_vcf,final_filtered_vcf,final_tsv,vep_summary,tumor_snv_bam_readcount_tsv,tumor_indel_bam_readcount_tsv,normal_snv_bam_readcount_tsv,normal_indel_bam_readcount_tsv,intervals_antitarget,intervals_target,normal_antitarget_coverage,normal_target_coverage,reference_coverage,cn_diagram,cn_scatter_plot,tumor_antitarget_coverage,tumor_target_coverage,tumor_bin_level_ratios,tumor_segmented_ratios,diploid_variants,somatic_variants,all_candidates,small_candidates,tumor_only_variants,somalier_concordance_metrics,somalier_concordance_statistics]
    germline:
        run: germline_exome_hla_typing.cwl
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
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: gatk_haplotypecaller_intervals
            vep_cache_dir: vep_cache_dir
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            vep_custom_annotations: vep_custom_annotations
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            qc_minimum_mapping_quality: qc_minimum_mapping_quality
            qc_minimum_base_quality: qc_minimum_base_quality
            optitype_name: optitype_name
        out:
            [cram,mark_duplicates_metrics,insert_size_metrics,insert_size_histogram,alignment_summary_metrics,hs_metrics,per_target_coverage_metrics,per_target_hs_metrics,per_base_coverage_metrics,per_base_hs_metrics,summary_hs_metrics,flagstats,verify_bam_id_metrics,verify_bam_id_depth,gvcf,final_vcf,filtered_vcf,vep_summary,optitype_tsv,optitype_plot]

    phase_vcf:
        run: ../subworkflows/phase_vcf.cwl
        in:
            somatic_vcf: somatic/final_vcf
            germline_vcf: germline/final_vcf
            reference: reference
            reference_dict: reference_dict
            bam: somatic/tumor_cram
            normal_sample_name: normal_sample_name
            tumor_sample_name: tumor_sample_name
        out:
            [phased_vcf]
    extract_alleles:
        run: ../tools/extract_hla_alleles.cwl
        in:
            allele_file: germline/optitype_tsv
        out:
            [allele_string]
    hla_consensus:
        run: ../tools/hla_consensus.cwl
        in:
            optitype_calls: extract_alleles/allele_string
            clinical_calls: clinical_calls
        out:
            [consensus_alleles, hla_call_files]
    pvacseq:
        run: ../subworkflows/pvacseq.cwl
        in:
            detect_variants_vcf: somatic/final_vcf
            sample_name: tumor_sample_name
            normal_sample_name: normal_sample_name
            rnaseq_bam: rnaseq/final_bam
            reference_fasta: reference
            readcount_minimum_base_quality: readcount_minimum_base_quality
            readcount_minimum_mapping_quality: readcount_minimum_mapping_quality
            gene_expression_file: rnaseq/gene_abundance
            transcript_expression_file: rnaseq/transcript_abundance_tsv
            alleles: hla_consensus/consensus_alleles
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
            variants_to_table_fields: variants_to_table_fields
            variants_to_table_genotype_fields: variants_to_table_genotype_fields
            vep_to_table_fields: vep_to_table_fields
        out:
            [annotated_vcf, annotated_tsv, mhc_i_all_epitopes, mhc_i_filtered_epitopes, mhc_i_ranked_epitopes, mhc_ii_all_epitopes, mhc_ii_filtered_epitopes, mhc_ii_ranked_epitopes, combined_all_epitopes, combined_filtered_epitopes, combined_ranked_epitopes]
