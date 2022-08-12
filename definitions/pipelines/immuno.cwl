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
        label: "tumor_sequence: MT sequencing data and readgroup information"
        doc: |
          tumor_sequence represents the sequencing data for the MT sample as either FASTQs or BAMs with
          accompanying readgroup information. The readgroup field should contain an entire read group header
          line, as described in the SAM file specification. This is a list of strings, beginning with @RG and
          followed by key:value pairs; each element of the list should be separated by a tab (\t). Keys ID and
          SM are required; see below for a formatting example:
          readgroup: "@RG\tID:xxx\tSM:xx"
          sequence:
            fastq1:
                class: File
                path: /path/to/reads1.fastq
            fastq2:
                class: File
                path: /path/to/reads2.fastq
            OR
            bam:
                class: File
                path: /path/to/reads.bam
    tumor_filename:
        type: string?
        default: 'tumor'
        label: "tumor/MT aligned bam filename"
        doc: |
          the filename to be used for bam files produced by the pipeline containing aligned
          tumor/mutant reads
    normal_sequence:
        type: ../types/sequence_data.yml#sequence_data[]
        label: "normal_sequence: WT sequencing data and readgroup information"
        doc: |
          normal_sequence represents the sequencing data for the WT sample as either FASTQs or BAMs with
          accompanying readgroup information. The readgroup field should contain an entire read group header
          line, as described in the SAM file specification. This is a list of strings, beginning with @RG and
          followed by key:value pairs; each element of the list should be separated by a tab (\t). Keys ID and
          SM are required; see below for a formatting example:
          readgroup: "@RG\tID:xxx\tSM:xx"
          sequence:
            fastq1:
                class: File
                path: /path/to/reads1.fastq
            fastq2:
                class: File
                path: /path/to/reads2.fastq
            OR
            bam:
                class: File
                path: /path/to/reads.bam
    normal_filename:
        type: string?
        default: 'normal'
        label: "normal/WT aligned bam filename"
        doc: |
          the filename to be used for bam files produced by the pipeline containing aligned
          normal/wild-type reads
    bqsr_known_sites:
        type: File[]
        secondaryFiles: [.tbi]
        label: "bqsr_known_sites: One or more databases of known polymorphic sites used to exclude regions around known polymorphisms from analysis."
        doc: |
          Known polymorphic indels recommended by GATK for a variety of
          tools including the BaseRecalibrator. This is part of the GATK resource
          bundle available at http://www.broadinstitute.org/gatk/guide/article?id=1213
          File should be in vcf format, and tabix indexed.
    bqsr_intervals:
        type: string[]
        label: "bqsr_intervals: Array of strings specifying regions for base quality score recalibration"
        doc: |
          bqsr_intervals provides an array of genomic intervals for which to apply
          GATK base quality score recalibrations. Typically intervals are given
          for the entire chromosome (chr1, chr2, etc.), these names should match
          the format in the reference file.
    bait_intervals:
        type: File
        label: "bait_intervals: interval_list file of baits used in the sequencing experiment"
        doc: |
          bait_intervals is an interval_list corresponding to the baits used in sequencing reagent.
          These are essentially coordinates for regions you were able to design probes for in the reagent.
          Typically the reagent provider has this information available in bed format and it can be
          converted to an interval_list with Picard BedToIntervalList. AstraZeneca also maintains a repo
          of baits for common sequencing reagents available at https://github.com/AstraZeneca-NGS/reference_data
    target_intervals:
        type: File
        label: "target_intervals: interval_list file of targets used in the sequencing experiment"
        doc: |
          target_intervals is an interval_list corresponding to the targets for the capture reagent.
          BED files with this information can be converted to interval_lists with Picard BedToIntervalList.
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
        doc: "Common mutations in cancer that will be genotyped and passed through into the merged VCF if they have even low-level evidence of a mutation (by default, marked with filter DOCM_ONLY)"
    filter_docm_variants:
        type: boolean?
        default: true
        doc: "Determines whether variants found only via genotyping of DOCM sites will be filtered (as DOCM_ONLY) or passed through as variant calls"
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
    validated_variants:
        type: File?
        secondaryFiles: [.tbi]
        doc: "An optional VCF with variants that will be flagged as 'VALIDATED' if found in this pipeline's main output VCF"

    #germline inputs
    gvcf_gq_bands:
        type: string[]
    gatk_haplotypecaller_intervals:
        type:
            type: array
            items:
                type: array
                items: string
    ploidy:
        type: int?
    optitype_name:
        type: string?

    #phase_vcf inputs
    reference_dict:
        type: File

    clinical_mhc_classI_alleles:
        type: string[]?
        label: "Clinical HLA typing results, limited to MHC Class I alleles; element format: HLA-X*01:02[/HLA-X...]"
        doc: "used to provide clinical HLA typing results in the format HLA-X*01:02[/HLA-X...] when available."
    clinical_mhc_classII_alleles:
        type: string[]?
        label: "Clinical HLA typing results, limited to MHC Class II alleles"
        doc: "used to provide clinical HLA typing results; separated from class I due to nomenclature inconsistencies"
    hla_source_mode:
        type:
            type: enum
            symbols: ["consensus", "clinical_only"]
        label: "Source for HLA types used for epitope prediction: in silico and clinical, or just clinical"
        doc: |
            Control whether HLA types passed to pvacseq should be a consensus of optitype predictions and clinical calls,
            if provided, or if only clinical calls should be used. In this case, optitype predictions and mismatches
            between optitype and clinical calls will still be reported. Selecting clinical_only without providing
            clinical calls will result in an error.

    #pvacseq inputs
    readcount_minimum_base_quality:
        type: int?
    readcount_minimum_mapping_quality:
        type: int?
    prediction_algorithms:
        type: string[]
    epitope_lengths_class_i:
        type: int[]?
    epitope_lengths_class_ii:
        type: int[]?
    binding_threshold:
        type: int?
    percentile_threshold:
        type: int?
    allele_specific_binding_thresholds:
        type: boolean?
    minimum_fold_change:
        type: float?
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
    expression_tool:
        type: string?
        default: 'kallisto'
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
        label: "net_chop_method: NetChop prediction method to use ('cterm' for C term 3.0, '20s' for 20S 3.0)"
        doc: |
           net_chop_method is used to specify which NetChop prediction method to use ("cterm" for C term 3.0, "20s" for 20S 3.0).
           C-term 3.0 is trained with publicly available MHC class I ligands and the authors believe that is performs best in predicting the
           boundaries of CTL epitopes. 20S is trained with in vitro degradation data.
    net_chop_threshold:
        type: float?
        label: "net_chop_threshold: NetChop prediction threshold"
        doc: |
          net_chop_threshold specifies the threshold to use for NetChop prediction; increasing the threshold results in better specificity, but worse sensitivity.
    netmhc_stab:
        type: boolean?
        label: "netmhc_stab: sets an option whether to run  NetMHCStabPan or not"
        doc: |
          netmhc_stab sets an option that decides whether it will run NetMHCStabPan after all filtering and add stability predictions to predicted epitopes.
    run_reference_proteome_similarity:
        type: boolean?
        label: "run_reference_proteome_similarity: sets an option whether to run reference proteome similarity or not"
        doc: |
          run_reference_proteome_similarity sets an option that decides whether it will run reference proteome similarity after all filtering and BLAST peptide sequences against the reference proteome to see if they appear elsewhere in the proteome.
    blastp_db:
        type:
            - "null"
            - type: enum
              symbols: ["refseq_select_prot", "refseq_protein"]
        label: "blastp_db: sets the reference proteome database to use with BLASTp"
        doc: |
          blastp_db sets the reference proteome database to use with BLASTp when enabling run_reference_proteome_similarity
    pvacseq_threads:
        type: int?
        label: "pvacseq_threads: Number of threads to use for parallelizing pvacseq prediction"
        doc: |
          pvacseq_threads specifies the number of threads to use for parallelizing peptide-MHC binding prediction calls.
    tumor_sample_name:
        type: string
        label: "tumor_sample_name: Name of the tumor sample"
        doc: |
          tumor_sample_name is the name of the tumor sample being processed. When processing a multi-sample VCF the sample name must be a sample ID in the input VCF #CHROM header line.
    normal_sample_name:
        type: string
        label: "normal_sample_name: Name of the normal sample"
        doc: |
          normal_sample_name is the name of the normal sample to use for phasing of germline variants.
    tumor_purity:
        type: float?

    #pvacfuse inputs
    iedb_retries:
        type: int?
    pvacfuse_keep_tmp_files:
        type: boolean?

    #FDA metrics inputs
    reference_genome_name:
        type: string?
    dna_sequencing_platform:
        type: string?
    dna_sequencing_instrument:
        type: string?
    dna_sequencing_kit:
        type: string?
    dna_sequencing_type:
        type: string?
    dna_single_or_paired_end:
        type: string?
    normal_dna_spike_in_error_rate:
        type: string?
    tumor_dna_spike_in_error_rate:
        type: string?
    normal_dna_total_DNA:
        type: string?
    tumor_dna_total_DNA:
        type: string?

outputs:
    tumor_cram:
        type: File
        outputSource: somatic/tumor_cram
        label: "Sorted CRAM from tumor DNA"
        doc: |
          Sorted CRAM file of sequencing read alignments by bwa-mem from a tumor DNA sample with duplicate reads tagged
    tumor_mark_duplicates_metrics:
        type: File
        outputSource: somatic/tumor_mark_duplicates_metrics
        label: "Sequencing duplicate metrics from tumor DNA"
        doc: |
          Duplication metrics on duplicate sequencing reads from a tumor DNA sample, identified by picard MarkDuplicates tool
    tumor_insert_size_metrics:
        type: File
        outputSource: somatic/tumor_insert_size_metrics
        label: "Paired-end sequencing diagnosis/quality metrics from tumor DNA"
        doc: |
          Diagnosis/quality metrics including the insert size distribution and read orientation of the paired-end libraries from a tumor DNA sample
    tumor_alignment_summary_metrics:
        type: File
        outputSource: somatic/tumor_alignment_summary_metrics
        label: "Sequencign alignment summary from tumor DNA"
        doc: |
          Diagnosis/quality metrics summarizing the quality of sequencing read alignments from a tumor DNA sample, reported by the picard CollectAlignmentSummaryMetrics tool
    tumor_hs_metrics:
        type: File
        outputSource: somatic/tumor_hs_metrics
        label: "Sequencing coverage summary of target intervals from tumor DNA"
        doc: |
          Diagnosis/quality metrics specific for sequencing data generated through hybrid-selection (e.g. whole exome) from a tumor DNA sample, for example to assess target coverage of WES
    tumor_per_target_coverage_metrics:
        type: File[]
        outputSource: somatic/tumor_per_target_coverage_metrics
        label: "Sequencing per-target coverage summary of target intervals from tumor DNA"
        doc: |
          Diagnosis/quality metrics showing detailed sequencing coverage per target interval (optional, 59 genes recommended by ACMG for clinical exome and genome sequencing for example) from a tumor DNA sample
    tumor_per_target_hs_metrics:
        type: File[]
        outputSource: somatic/tumor_per_target_hs_metrics
        label: "Sequencing coverage summary of target intervals from tumor DNA"
        doc: |
          Diagnosis/quality metrics for sequencing coverage for target intervals (optional, 59 genes recommended by ACMG for clinical exome and genome sequencing for example) from a tumor DNA sample
    tumor_per_base_coverage_metrics:
        type: File[]
        outputSource: somatic/tumor_per_base_coverage_metrics
        label: "Sequencing per-base coverage summary at target sites from tumor DNA"
        doc: |
          Diagnosis/quality metrics showing detailed sequencing coverage per target site (optional, known variant sites of clinical significance from ClinVar for example) from a tumor DNA sample
    tumor_per_base_hs_metrics:
        type: File[]
        outputSource: somatic/tumor_per_base_hs_metrics
        label: "Sequencing coverage summary at target sites from tumor DNA"
        doc: |
          Diagnosis/quality metrics for sequencing coverage at target sites (optional, known variant sites of clinical significance from ClinVar for example) from a tumor DNA sample
    tumor_summary_hs_metrics:
        type: File[]
        outputSource: somatic/tumor_summary_hs_metrics
    tumor_flagstats:
        type: File
        outputSource: somatic/tumor_flagstats
        label: "Sequencing count metrics based on SAM FLAG field from tumor sample"
        doc: |
          Summary with the count numbers of alignments for each FLAG type from a tumor DNA sample, including 13 categories based on the bit flags in the FLAG field
    tumor_verify_bam_id_metrics:
        type: File
        outputSource: somatic/tumor_verify_bam_id_metrics
        label: "Sequencing quality assessment metric for tumor sample contamination"
        doc: |
          verifyBamID output files containing the contamination estimate in a tumor DNA sample, across all readGroups and per readGroup separately
    tumor_verify_bam_id_depth:
        type: File
        outputSource: somatic/tumor_verify_bam_id_depth
        label: "Sequencing quality assessment metric for tumor sample genotyping"
        doc: |
          verifyBamID output files showing the sequencing depth distribution at the marker positions from Omni genotype data with a tumor DNA sample, across all readGroups and per readGroup separately
    normal_cram:
        type: File
        outputSource: somatic/normal_cram
        label: "Sorted CRAM from normal DNA"
        doc: |
          Sorted CRAM file of sequencing read alignments by bwa-mem from a normal DNA sample with duplicate reads tagged
    normal_mark_duplicates_metrics:
        type: File
        outputSource: somatic/normal_mark_duplicates_metrics
        label: "Sequencing duplicate metrics from normal DNA"
        doc: |
          Duplication metrics on duplicate sequencing reads from a normal DNA sample, identified by picard MarkDuplicates tool
    normal_insert_size_metrics:
        type: File
        outputSource: somatic/normal_insert_size_metrics
        label: "Paired-end sequencing diagnosis/quality metrics from normal DNA"
        doc: |
          Diagnosis/quality metrics including the insert size distribution and read orientation of the paired-end libraries from a normal DNA sample
    normal_alignment_summary_metrics:
        type: File
        outputSource: somatic/normal_alignment_summary_metrics
        label: "Sequencign alignment summary from normal DNA"
        doc: |
          Diagnosis/quality metrics summarizing the quality of sequencing read alignments from a normal DNA sample, reported by the picard CollectAlignmentSummaryMetrics tool
    normal_hs_metrics:
        type: File
        outputSource: somatic/normal_hs_metrics
        label: "Sequencing coverage summary of target intervals from normal DNA"
        doc: |
          Diagnosis/quality metrics specific for sequencing data generated through hybrid-selection (e.g. whole exome) from a normal DNA sample, for example to assess target coverage
    normal_per_target_coverage_metrics:
        type: File[]
        outputSource: somatic/normal_per_target_coverage_metrics
        label: "Sequencing per-target coverage summary of target intervals from normal DNA"
        doc: |
          Diagnosis/quality metrics showing detailed sequencing coverage per target interval (optional, 59 genes recommended by ACMG for clinical exome and genome sequencing for example) from a normal DNA sample
    normal_per_target_hs_metrics:
        type: File[]
        outputSource: somatic/normal_per_target_hs_metrics
        label: "Sequencing coverage summary of target intervals from normal DNA"
        doc: |
          Diagnosis/quality metrics for sequencing coverage for target intervals (optional, 59 genes recommended by ACMG for clinical exome and genome sequencing for example) from a normal DNA sample
    normal_per_base_coverage_metrics:
        type: File[]
        outputSource: somatic/normal_per_base_coverage_metrics
        label: "Sequencing per-base coverage summary at target sites from normal DNA"
        doc: |
          Diagnosis/quality metrics showing detailed sequencing coverage per target site (optional, known variant sites of clinical significance from ClinVar for example) from a normal DNA sample
    normal_per_base_hs_metrics:
        type: File[]
        outputSource: somatic/normal_per_base_hs_metrics
        label: "Sequencing coverage summary at target sites from normal DNA"
        doc: |
          Diagnosis/quality metrics for sequencing coverage at target sites (optional, known variant sites of clinical significance from ClinVar for example) from a normal DNA sample
    normal_summary_hs_metrics:
        type: File[]
        outputSource: somatic/normal_summary_hs_metrics
    normal_flagstats:
        type: File
        outputSource: somatic/normal_flagstats
        label: "Sequencing count metrics based on SAM FLAG field from normal sample"
        doc: |
          Summary with the count numbers of alignments for each FLAG type from a normal DNA sample, including 13 categories based on the bit flags in the FLAG field
    normal_verify_bam_id_metrics:
        type: File
        outputSource: somatic/normal_verify_bam_id_metrics
        label: "Sequencing quality assessment metric for normal sample contamination"
        doc: |
          verifyBamID output files containing the contamination estimate in a normal DNA sample, across all readGroups and per readGroup separately
    normal_verify_bam_id_depth:
        type: File
        outputSource: somatic/normal_verify_bam_id_depth
        label: "Sequencing quality assessment metric for normal sample genotyping"
        doc: |
          verifyBamID output files showing the sequencing depth distribution at the marker positions from Omni genotype data with a normal DNA sample, across all readGroups and per readGroup separately
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
    germline_raw_vcf:
        type: File
        outputSource: germline/raw_vcf
        secondaryFiles: [.tbi]
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
    pvacseq_predictions:
        type: Directory
        outputSource: pvacseq/pvacseq_predictions
    pvacfuse_predictions:
        type: Directory
        outputSource: pvacfuse/pvacfuse_predictions

    unaligned_normal_dna_fastqc_data:
        type: File[]
        outputSource: fda_metrics/unaligned_normal_dna_fastqc_data
    unaligned_normal_dna_table_metrics:
        type: File
        outputSource: fda_metrics/unaligned_normal_dna_table_metrics
    unaligned_normal_dna_md5sums:
        type: File
        outputSource: fda_metrics/unaligned_normal_dna_md5sums
    unaligned_normal_dna_table1:
        type: File
        outputSource: fda_metrics/unaligned_normal_dna_table1

    unaligned_tumor_dna_fastqc_data:
        type: File[]
        outputSource: fda_metrics/unaligned_tumor_dna_fastqc_data
    unaligned_tumor_dna_table_metrics:
        type: File
        outputSource: fda_metrics/unaligned_tumor_dna_table_metrics
    unaligned_tumor_dna_md5sums:
        type: File
        outputSource: fda_metrics/unaligned_tumor_dna_md5sums
    unaligned_tumor_dna_table1:
        type: File
        outputSource: fda_metrics/unaligned_tumor_dna_table1

    aligned_normal_dna_fastqc_data:
        type: File[]
        outputSource: fda_metrics/aligned_normal_dna_fastqc_data
    aligned_normal_dna_table_metrics:
        type: File
        outputSource: fda_metrics/aligned_normal_dna_table_metrics
    aligned_normal_dna_md5sums:
        type: File
        outputSource: fda_metrics/aligned_normal_dna_md5sums
    aligned_normal_dna_table2:
        type: File
        outputSource: fda_metrics/aligned_normal_dna_table2

    aligned_tumor_dna_fastqc_data:
        type: File[]
        outputSource: fda_metrics/aligned_tumor_dna_fastqc_data
    aligned_tumor_dna_table_metrics:
        type: File
        outputSource: fda_metrics/aligned_tumor_dna_table_metrics
    aligned_tumor_dna_md5sums:
        type: File
        outputSource: fda_metrics/aligned_tumor_dna_md5sums
    aligned_tumor_dna_table2:
        type: File
        outputSource: fda_metrics/aligned_tumor_dna_table2

steps:
    somatic:
        run: somatic_exome.cwl
        in:
            reference: reference
            tumor_sequence: tumor_sequence
            tumor_name: tumor_filename
            normal_sequence: normal_sequence
            normal_name: normal_filename
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
            qc_minimum_mapping_quality: qc_minimum_mapping_quality
            qc_minimum_base_quality: qc_minimum_base_quality
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
            validated_variants: validated_variants
        out:
            [tumor_cram,tumor_mark_duplicates_metrics,tumor_insert_size_metrics,tumor_alignment_summary_metrics,tumor_hs_metrics,tumor_per_target_coverage_metrics,tumor_per_target_hs_metrics,tumor_per_base_coverage_metrics,tumor_per_base_hs_metrics,tumor_summary_hs_metrics,tumor_flagstats,tumor_verify_bam_id_metrics,tumor_verify_bam_id_depth,normal_cram,normal_mark_duplicates_metrics,normal_insert_size_metrics,normal_alignment_summary_metrics,normal_hs_metrics,normal_per_target_coverage_metrics,normal_per_target_hs_metrics,normal_per_base_coverage_metrics,normal_per_base_hs_metrics,normal_summary_hs_metrics,normal_flagstats,normal_verify_bam_id_metrics,normal_verify_bam_id_depth,mutect_unfiltered_vcf,mutect_filtered_vcf,strelka_unfiltered_vcf,strelka_filtered_vcf,varscan_unfiltered_vcf,varscan_filtered_vcf,pindel_unfiltered_vcf,pindel_filtered_vcf,docm_filtered_vcf,final_vcf,final_filtered_vcf,final_tsv,vep_summary,tumor_snv_bam_readcount_tsv,tumor_indel_bam_readcount_tsv,normal_snv_bam_readcount_tsv,normal_indel_bam_readcount_tsv,intervals_antitarget,intervals_target,normal_antitarget_coverage,normal_target_coverage,reference_coverage,cn_diagram,cn_scatter_plot,tumor_antitarget_coverage,tumor_target_coverage,tumor_bin_level_ratios,tumor_segmented_ratios,diploid_variants,somatic_variants,all_candidates,small_candidates,tumor_only_variants,somalier_concordance_metrics,somalier_concordance_statistics]
    germline:
        run: germline_exome_hla_typing.cwl
        in:
            reference: reference
            sequence: normal_sequence
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
            gvcf_gq_bands: gvcf_gq_bands
            intervals: gatk_haplotypecaller_intervals
            ploidy: ploidy
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
            [cram,mark_duplicates_metrics,insert_size_metrics,insert_size_histogram,alignment_summary_metrics,hs_metrics,per_target_coverage_metrics,per_target_hs_metrics,per_base_coverage_metrics,per_base_hs_metrics,summary_hs_metrics,flagstats,verify_bam_id_metrics,verify_bam_id_depth,raw_vcf,final_vcf,filtered_vcf,vep_summary,optitype_tsv,optitype_plot]

    fda_metrics:
        run: ../subworkflows/generate_fda_metrics.cwl
        in:
            reference: reference
            unaligned_normal_dna: normal_sequence
            unaligned_tumor_dna: tumor_sequence
            aligned_normal_dna: somatic/normal_cram
            aligned_tumor_dna: somatic/tumor_cram
            normal_alignment_summary_metrics: somatic/normal_alignment_summary_metrics
            normal_duplication_metrics: somatic/normal_mark_duplicates_metrics
            normal_insert_size_metrics: somatic/normal_insert_size_metrics
            normal_hs_metrics: somatic/normal_hs_metrics
            normal_flagstat: somatic/normal_flagstats
            tumor_alignment_summary_metrics: somatic/tumor_alignment_summary_metrics
            tumor_duplication_metrics: somatic/tumor_mark_duplicates_metrics
            tumor_insert_size_metrics: somatic/tumor_insert_size_metrics
            tumor_hs_metrics: somatic/tumor_hs_metrics
            tumor_flagstat: somatic/tumor_flagstats
            reference_genome: reference_genome_name
            dna_sequencing_platform: dna_sequencing_platform
            dna_sequencing_instrument: dna_sequencing_instrument
            dna_sequencing_kit: dna_sequencing_kit
            dna_sequencing_type: dna_sequencing_type
            dna_single_or_paired_end: dna_single_or_paired_end
            normal_dna_spike_in_error_rate: normal_dna_spike_in_error_rate
            tumor_dna_spike_in_error_rate: tumor_dna_spike_in_error_rate
            normal_dna_total_DNA: normal_dna_total_DNA
            tumor_dna_total_DNA: tumor_dna_total_DNA
            normal_dna_sample_name: normal_sample_name
            tumor_dna_sample_name: tumor_sample_name
        out:
            [unaligned_normal_dna_fastqc_data,unaligned_normal_dna_table_metrics,unaligned_normal_dna_md5sums,unaligned_normal_dna_table1,unaligned_tumor_dna_fastqc_data,unaligned_tumor_dna_table_metrics,unaligned_tumor_dna_md5sums,unaligned_tumor_dna_table1,aligned_normal_dna_fastqc_data,aligned_normal_dna_table_metrics,aligned_normal_dna_md5sums,aligned_normal_dna_table2,aligned_tumor_dna_fastqc_data,aligned_tumor_dna_table_metrics,aligned_tumor_dna_md5sums,aligned_tumor_dna_table2]

    phase_vcf:
        run: ../subworkflows/phase_vcf.cwl
        in:
            somatic_vcf: somatic/final_filtered_vcf
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
            hla_source_mode: hla_source_mode
            optitype_hla_alleles: extract_alleles/allele_string
            clinical_mhc_classI_alleles: clinical_mhc_classI_alleles
            clinical_mhc_classII_alleles: clinical_mhc_classII_alleles
        out:
            [consensus_alleles, hla_call_files]
    intersect_passing_variants:
        run: ../tools/intersect_known_variants.cwl
        in:
            vcf: somatic/final_filtered_vcf
            validated_variants: validated_variants
        out:
            [validated_and_pipeline_vcf]
    pvacseq:
        run: ../subworkflows/pvacseq.cwl
        in:
            detect_variants_vcf: intersect_passing_variants/validated_and_pipeline_vcf
            sample_name: tumor_sample_name
            normal_sample_name: normal_sample_name
            reference_fasta: reference
            readcount_minimum_base_quality: readcount_minimum_base_quality
            readcount_minimum_mapping_quality: readcount_minimum_mapping_quality
            alleles: hla_consensus/consensus_alleles
            prediction_algorithms: prediction_algorithms
            epitope_lengths_class_i: epitope_lengths_class_i
            epitope_lengths_class_ii: epitope_lengths_class_ii
            binding_threshold: binding_threshold
            percentile_threshold: percentile_threshold
            allele_specific_binding_thresholds: allele_specific_binding_thresholds
            minimum_fold_change: minimum_fold_change
            top_score_metric: top_score_metric
            additional_report_columns: additional_report_columns
            expression_tool: expression_tool
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
            run_reference_proteome_similarity: run_reference_proteome_similarity
            blastp_db: blastp_db
            n_threads: pvacseq_threads
            variants_to_table_fields: variants_to_table_fields
            variants_to_table_genotype_fields: variants_to_table_genotype_fields
            vep_to_table_fields: vep_to_table_fields
            tumor_purity: tumor_purity
        out:
            [annotated_vcf, annotated_tsv, pvacseq_predictions]
