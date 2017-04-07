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
        type: File
    synonyms_file:
        type: File?
outputs:
    final_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [.tbi]
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
    combine:
        run: combine_variants.cwl
        in:
            reference: reference
            mutect_vcf: mutect/merged_vcf
            strelka_vcf: strelka/merged_vcf
            varscan_vcf: varscan/merged_vcf
            pindel_vcf: pindel/merged_vcf
        out:
            [combined_vcf]
    filter:
        run: ../fp_filter/workflow.cwl
        in:
            reference: reference
            cram: tumor_cram
            vcf: combine/combined_vcf
        out:
            [filtered_vcf]
    fp_bgzip:
        run: bgzip.cwl
        in:
            file: filter/filtered_vcf
        out:
            [bgzipped_file]
    fp_index:
        run: index.cwl
        in:
            vcf: fp_bgzip/bgzipped_file
        out:
            [indexed_vcf]
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
    combine_docm:
        run: combine_docm.cwl
        in: 
            reference: reference
            filter_vcf: fp_index/indexed_vcf
            docm_vcf: docm/merged_vcf
        out:
            [combine_docm_vcf]
    annotate_variants:
        run: vep.cwl
        in:
            vcf: combine_docm/combine_docm_vcf
            cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
        out:
            [annotated_vcf]
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
