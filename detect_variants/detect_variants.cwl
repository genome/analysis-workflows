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
    tumor_bam:
        type: File
        secondaryFiles: [".bai"]
    normal_bam:
        type: File
        secondaryFiles: [".bai"]
    interval_list:
        type: File
    dbsnp_vcf:
        type: File?
        secondaryFiles: [".tbi"]
    cosmic_vcf:
        type: File?
        secondaryFiles: [".tbi"]
    strelka_exome_mode:
        type: boolean
    pindel_chromosome_list:
        type:
            type: array
            items: string
    pindel_insert_size:
        type: int
        default: 400
outputs:
    final_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: [".tbi"]
steps:
    mutect:
        run: ../mutect/workflow.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: interval_list
            dbsnp_vcf: dbsnp_vcf
            cosmic_vcf: cosmic_vcf
        out:
            [merged_vcf]
    strelka:
        run: ../strelka/workflow.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: interval_list
            exome_mode: strelka_exome_mode
        out:
            [merged_vcf]
    varscan:
        run: ../varscan/workflow.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: interval_list
        out:
            [merged_vcf]
    pindel:
        run: ../pindel/workflow.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: interval_list
            chromosome_list: pindel_chromosome_list
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
        run: fp_filter.cwl
        in:
            reference: reference
            bam: tumor_bam
            vcf: combine/combined_vcf
        out:
            [filtered_vcf]
    bgzip:
        run: bgzip.cwl
        in:
            file: filter/filtered_vcf
        out:
            [bgzipped_file]
    index:
        run: index.cwl
        in:
            vcf: bgzip/bgzipped_file
        out:
            [indexed_vcf]
