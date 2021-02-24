#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Detect Variants workflow"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    tumor_bam:
        type: File
        secondaryFiles: [.bai,^.bai]
    normal_bam:
        type: File
        secondaryFiles: [.bai,^.bai]
    roi_intervals:
        type: File
        label: "roi_intervals: regions of interest in which variants will be called"
        doc: "roi_intervals is a list of regions (in interval_list format) within which to call somatic variants"
    strelka_exome_mode:
        type: boolean
    strelka_cpu_reserved:
        type: int?
        default: 8
    readcount_minimum_base_quality:
        type: int?
    readcount_minimum_mapping_quality:
        type: int?
    scatter_count:
        type: int
        doc: "scatters each supported variant detector (varscan, pindel, mutect) into this many parallel jobs"
    varscan_strand_filter:
        type: int?
        default: 0
    varscan_min_coverage:
        type: int?
        default: 8
    varscan_min_var_freq:
        type: float?
        default: 0.1
    varscan_p_value:
        type: float?
        default: 0.99
    varscan_max_normal_freq:
        type: float?
    pindel_insert_size:
        type: int
        default: 400
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
    vep_plugins:
        type: string[]
        default: [Frameshift, Wildtype]
    filter_mapq0_threshold:
        type: float
        default: 0.15
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
    filter_minimum_depth:
        type: int
        default: 1
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
        default: []
    tumor_sample_name:
        type: string
    normal_sample_name:
        type: string
outputs:
    mutect_unfiltered_vcf:
        type: File
        outputSource: mutect/unfiltered_vcf
        secondaryFiles: [.tbi]
    mutect_filtered_vcf:
        type: File
        outputSource: mutect/filtered_vcf
        secondaryFiles: [.tbi]
    strelka_unfiltered_vcf:
        type: File
        outputSource: strelka/unfiltered_vcf
        secondaryFiles: [.tbi]
    strelka_filtered_vcf:
        type: File
        outputSource: strelka/filtered_vcf
        secondaryFiles: [.tbi]
    varscan_unfiltered_vcf:
        type: File
        outputSource: varscan/unfiltered_vcf
        secondaryFiles: [.tbi]
    varscan_filtered_vcf:
        type: File
        outputSource: varscan/filtered_vcf
        secondaryFiles: [.tbi]
    pindel_unfiltered_vcf:
        type: File
        outputSource: pindel/unfiltered_vcf
        secondaryFiles: [.tbi]
    pindel_filtered_vcf:
        type: File
        outputSource: pindel/filtered_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [.tbi]
    final_filtered_vcf:
        type: File
        outputSource: annotated_filter_index/indexed_vcf
        secondaryFiles: [.tbi]
    final_tsv:
        type: File
        outputSource: add_vep_fields_to_table/annotated_variants_tsv
    vep_summary:
        type: File
        outputSource: annotate_variants/vep_summary
    tumor_snv_bam_readcount_tsv:
        type: File
        outputSource: tumor_bam_readcount/snv_bam_readcount_tsv
    tumor_indel_bam_readcount_tsv:
        type: File
        outputSource: tumor_bam_readcount/indel_bam_readcount_tsv
    normal_snv_bam_readcount_tsv:
        type: File
        outputSource: normal_bam_readcount/snv_bam_readcount_tsv
    normal_indel_bam_readcount_tsv:
        type: File
        outputSource: normal_bam_readcount/indel_bam_readcount_tsv
steps:
    mutect:
        run: ../subworkflows/mutect.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: roi_intervals
            scatter_count: scatter_count
            tumor_sample_name: tumor_sample_name
        out:
            [unfiltered_vcf, filtered_vcf]
    strelka:
        run: ../subworkflows/strelka_and_post_processing.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: roi_intervals
            exome_mode: strelka_exome_mode
            cpu_reserved: strelka_cpu_reserved
            tumor_sample_name: tumor_sample_name
            normal_sample_name: normal_sample_name
        out:
            [unfiltered_vcf, filtered_vcf]
    varscan:
        run: ../subworkflows/varscan_pre_and_post_processing.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: roi_intervals
            scatter_count: scatter_count
            strand_filter: varscan_strand_filter
            min_coverage: varscan_min_coverage
            min_var_freq: varscan_min_var_freq
            p_value: varscan_p_value
            max_normal_freq: varscan_max_normal_freq
            tumor_sample_name: tumor_sample_name
            normal_sample_name: normal_sample_name
        out:
            [unfiltered_vcf, filtered_vcf]
    pindel:
        run: ../subworkflows/pindel.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: roi_intervals
            scatter_count: scatter_count
            insert_size: pindel_insert_size
            tumor_sample_name: tumor_sample_name
            normal_sample_name: normal_sample_name
        out:
            [unfiltered_vcf, filtered_vcf]
    combine:
        run: ../tools/combine_variants.cwl
        in:
            reference: reference
            mutect_vcf: mutect/filtered_vcf
            strelka_vcf: strelka/filtered_vcf
            varscan_vcf: varscan/filtered_vcf
            pindel_vcf: pindel/filtered_vcf
        out:
            [combined_vcf]
    decompose:
        run: ../tools/vt_decompose.cwl
        in:
            vcf: combine/combined_vcf
        out:
            [decomposed_vcf]
    decompose_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: decompose/decomposed_vcf
        out:
            [indexed_vcf]
    annotate_variants:
        run: ../tools/vep.cwl
        in:
            vcf: decompose_index/indexed_vcf
            cache_dir: vep_cache_dir
            ensembl_assembly: vep_ensembl_assembly
            ensembl_version: vep_ensembl_version
            ensembl_species: vep_ensembl_species
            synonyms_file: synonyms_file
            coding_only: annotate_coding_only
            reference: reference
            pick: vep_pick
            plugins: vep_plugins
        out:
            [annotated_vcf, vep_summary]
    tumor_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: annotate_variants/annotated_vcf
            sample: tumor_sample_name
            reference_fasta: reference
            bam: tumor_bam
            min_base_quality: readcount_minimum_base_quality
            min_mapping_quality: readcount_minimum_mapping_quality
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    normal_bam_readcount:
        run: ../tools/bam_readcount.cwl
        in:
            vcf: annotate_variants/annotated_vcf
            sample: normal_sample_name
            reference_fasta: reference
            bam: normal_bam
            min_base_quality: readcount_minimum_base_quality
            min_mapping_quality: readcount_minimum_mapping_quality
        out:
            [snv_bam_readcount_tsv, indel_bam_readcount_tsv]
    add_tumor_bam_readcount_to_vcf:
        run: ../subworkflows/vcf_readcount_annotator.cwl
        in:
            vcf: annotate_variants/annotated_vcf
            snv_bam_readcount_tsv: tumor_bam_readcount/snv_bam_readcount_tsv
            indel_bam_readcount_tsv: tumor_bam_readcount/indel_bam_readcount_tsv
            data_type:
                default: 'DNA'
            sample_name: tumor_sample_name
        out:
            [annotated_bam_readcount_vcf]
    add_normal_bam_readcount_to_vcf:
        run: ../subworkflows/vcf_readcount_annotator.cwl
        in:
            vcf: add_tumor_bam_readcount_to_vcf/annotated_bam_readcount_vcf
            snv_bam_readcount_tsv: normal_bam_readcount/snv_bam_readcount_tsv
            indel_bam_readcount_tsv: normal_bam_readcount/indel_bam_readcount_tsv
            data_type:
                default: 'DNA'
            sample_name: normal_sample_name
        out:
            [annotated_bam_readcount_vcf]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: add_normal_bam_readcount_to_vcf/annotated_bam_readcount_vcf
        out:
            [indexed_vcf]
    filter_vcf:
        run: ../subworkflows/filter_vcf_nonhuman.cwl
        in: 
            vcf: index/indexed_vcf
            filter_mapq0_threshold: filter_mapq0_threshold
            filter_somatic_llr_threshold: filter_somatic_llr_threshold
            filter_somatic_llr_tumor_purity: filter_somatic_llr_tumor_purity
            filter_somatic_llr_normal_contamination_rate: filter_somatic_llr_normal_contamination_rate
            filter_minimum_depth: filter_minimum_depth
            tumor_bam: tumor_bam
            do_cle_vcf_filter: cle_vcf_filter
            reference: reference
            tumor_sample_name: tumor_sample_name
            normal_sample_name: normal_sample_name
        out: 
            [filtered_vcf]
    annotated_filter_bgzip:
        run: ../tools/bgzip.cwl
        in:
            file: filter_vcf/filtered_vcf
        out:
            [bgzipped_file]
    annotated_filter_index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: annotated_filter_bgzip/bgzipped_file
        out:
            [indexed_vcf]
    variants_to_table:
        run: ../tools/variants_to_table.cwl
        in:
            reference: reference
            vcf: annotated_filter_index/indexed_vcf
            fields: variants_to_table_fields
            genotype_fields: variants_to_table_genotype_fields
        out:
            [variants_tsv]
    add_vep_fields_to_table:
        run: ../tools/add_vep_fields_to_table.cwl
        in:
            vcf: annotated_filter_index/indexed_vcf
            vep_fields: vep_to_table_fields
            tsv: variants_to_table/variants_tsv
        out: [annotated_variants_tsv]
