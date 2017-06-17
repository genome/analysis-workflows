#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Detect Variants workflow"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type: File
        secondaryFiles: [".fai", "^.dict"]
    tumor_cram:
        type: File
        secondaryFiles: [^.crai,.crai]
    normal_cram:
        type: File
        secondaryFiles: [^.crai,.crai]
    interval_list:
        type: File
    dbsnp_vcf:
        type: File?
        secondaryFiles: [.tbi]
    cosmic_vcf:
        type: File?
        secondaryFiles: [.tbi]
    panel_of_normals_vcf:
        type: File?
        secondaryFiles: [.tbi]
    strelka_exome_mode:
        type: boolean
    mutect_scatter_count:
        type: int?
    mutect_artifact_detection_mode:
        type: boolean?
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
        default: 0.1
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
    coding_only:
        type: boolean?
    hgvs_annotation:
        type: boolean?
    variants_to_table_fields:
        type: string[]?
        default: [CHROM,POS,ID,REF,ALT,set,AC,AF]
    variants_to_table_genotype_fields:
        type: string[]?
        default: [GT,AD]
    vep_to_table_fields:
        type: string[]?
        default: [HGVSc,HGVSp]
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
    docm_unfiltered_vcf:
        type: File
        outputSource: docm/unfiltered_vcf
    docm_filtered_vcf:
        type: File
        outputSource: docm/filtered_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [.tbi]
    final_tsv:
        type: File
        outputSource: add_vep_fields_to_table/annotated_variants_tsv
    vep_summary:
        type: File
        outputSource: annotate_variants/vep_summary
    snv_bam_readcount:
        type: File
        outputSource: bam_readcount/snv_bam_readcount
    indel_bam_readcount:
        type: File
        outputSource: bam_readcount/indel_bam_readcount
steps:
    mutect:
        run: ../mutect/workflow.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            interval_list: interval_list
            dbsnp_vcf: dbsnp_vcf
            cosmic_vcf: cosmic_vcf
            max_alt_allele_in_normal_fraction: mutect_max_alt_allele_in_normal_fraction
            max_alt_alleles_in_normal_count: mutect_max_alt_alleles_in_normal_count
            scatter_count: mutect_scatter_count
            artifact_detection_mode: mutect_artifact_detection_mode
            panel_of_normals_vcf: panel_of_normals_vcf
        out:
            [unfiltered_vcf, filtered_vcf]
    strelka:
        run: ../strelka/workflow.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            interval_list: interval_list
            exome_mode: strelka_exome_mode
        out:
            [unfiltered_vcf, filtered_vcf]
    varscan:
        run: ../varscan/workflow.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            interval_list: interval_list
            strand_filter: varscan_strand_filter
            min_coverage: varscan_min_coverage
            min_var_freq: varscan_min_var_freq
            p_value: varscan_p_value
            max_normal_freq: varscan_max_normal_freq
        out:
            [unfiltered_vcf, filtered_vcf]
    pindel:
        run: ../pindel/workflow.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            interval_list: interval_list
            insert_size: pindel_insert_size
        out:
            [unfiltered_vcf, filtered_vcf]
    docm:
        run: ../docm/workflow.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            docm_vcf: docm_vcf
            interval_list: interval_list
        out:
            [unfiltered_vcf, filtered_vcf]
    combine:
        run: combine_variants.cwl
        in:
            reference: reference
            mutect_vcf: mutect/filtered_vcf
            strelka_vcf: strelka/filtered_vcf
            varscan_vcf: varscan/filtered_vcf
            pindel_vcf: pindel/filtered_vcf
            docm_vcf: docm/filtered_vcf
        out:
            [combined_vcf]
    annotate_variants:
        run: vep.cwl
        in:
            vcf: combine/combined_vcf
            cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            coding_only: coding_only
            hgvs: hgvs_annotation
            reference: reference
        out:
            [annotated_vcf, vep_summary]
    cram_to_bam:
        run: ../cram_to_bam/workflow.cwl
        in:
            cram: tumor_cram
            reference: reference
        out:
            [bam]
    bam_readcount:
        run: ../pvacseq/bam_readcount.cwl
        in:
            vcf: combine/combined_vcf
            sample:
                default: 'TUMOR'
            reference_fasta: reference
            bam: cram_to_bam/bam
        out:
            [snv_bam_readcount, indel_bam_readcount]
    bgzip:
        run: bgzip.cwl
        in:
            file: annotate_variants/annotated_vcf
        out:
            [bgzipped_file]
    index:
        run: index.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
    variants_to_table:
        run: variants_to_table.cwl
        in:
            reference: reference
            vcf: index/indexed_vcf
            fields: variants_to_table_fields
            genotype_fields: variants_to_table_genotype_fields
        out:
            [variants_tsv]
    add_vep_fields_to_table:
        run: add_vep_fields_to_table.cwl
        in:
            vcf: index/indexed_vcf
            vep_fields: vep_to_table_fields
            tsv: variants_to_table/variants_tsv
        out: [annotated_variants_tsv]
