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
    strelka_exome_mode:
        type: boolean
    mutect_scatter_count:
        type: int?
    mutect_artifact_detection_mode:
        type: boolean?
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
    hard_filter_vcf:
        type: boolean?
        default: true
    variants_to_table_fields:
        type: string[]?
        default: [CHROM,POS,ID,REF,ALT,set,AC,AF]
    variants_to_table_genotype_fields:
        type: string[]?
        default: [GT,AD]
outputs:
    final_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [.tbi]
    final_tsv:
        type: File
        outputSource: variants_to_table/variants_tsv
    vep_summary:
        type: File
        outputSource: annotate_variants/vep_summary
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
            scatter_count: mutect_scatter_count
            artifact_detection_mode: mutect_artifact_detection_mode
        out:
            [merged_vcf]
    strelka:
        run: ../strelka/workflow.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            interval_list: interval_list
            exome_mode: strelka_exome_mode
        out:
            [merged_vcf]
    varscan:
        run: ../varscan/workflow.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            interval_list: interval_list
        out:
            [merged_vcf]
    pindel:
        run: ../pindel/workflow.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            interval_list: interval_list
            insert_size: pindel_insert_size
        out:
            [merged_vcf]
    docm:
        run: ../docm/workflow.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            docm_vcf: docm_vcf
            interval_list: interval_list
        out:
            [merged_vcf]
    combine:
        run: combine_variants.cwl
        in:
            reference: reference
            mutect_vcf: mutect/merged_vcf
            strelka_vcf: strelka/merged_vcf
            varscan_vcf: varscan/merged_vcf
            pindel_vcf: pindel/merged_vcf
            docm_vcf: docm/merged_vcf
        out:
            [combined_vcf]
    hard_filter:
        run: select_variants.cwl
        in:
            reference: reference
            vcf: combine/combined_vcf
            exclude_filtered: hard_filter_vcf
        out:
            [filtered_vcf]
    annotate_variants:
        run: vep.cwl
        in:
            vcf: hard_filter/filtered_vcf
            cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            coding_only: coding_only
        out:
            [annotated_vcf, vep_summary]
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
