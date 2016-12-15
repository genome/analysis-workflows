#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "somatic pipeline"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type: File
        secondaryFiles: [.fai, .bwt, .sa, .ann, .amb, .pac, ^.dict]
    normal_bams:
         type: File[]
    normal_readgroups:
          type: string[]
    tumor_bams:
         type: File[]
    tumor_readgroups:
          type: string[]
    alignment_dbsnp:
        type: File
        secondaryFiles: [.tbi]
    mills:
        type: File
        secondaryFiles: [.tbi]
    known_indels:
        type: File
        secondaryFiles: [.tbi]
    interval_list:
        type: File
    strelka_config:
        type: File
    dv_dbsnp_vcf:
        type: File?
        secondaryFiles: .idx
    cosmic_vcf:
        type: File?
        secondaryFiles: .idx
outputs:
    final_vcf:
        type: File
        outputSource: detect_variants/final_vcf
        secondaryFile: .tbi
steps:
    align_normal:
        run: unaligned_bam_to_bqsr_and_vcf/workflow.cwl
        in:
            bams: normal_bams
            readgroups: normal_readgroups
            reference: reference
            dbsnp: alignment_dbsnp
            mills: mills
            known_indels: known_indels
        out:
            [final_bam]
    align_tumor:
        run: unaligned_bam_to_bqsr_and_vcf/workflow.cwl
        in:
            bams: tumor_bams
            readgroups: tumor_readgroups
            reference: reference
            dbsnp: alignment_dbsnp
            mills: mills
            known_indels: known_indels
        out:
            [final_bam]
    detect_variants:
        run: detect_variants/detect_variants.cwl
        in:
            reference: reference
            tumor_bam: align_tumor/final_bam
            normal_bam: align_normal/final_bam
            interval_list: interval_list
            strelka_config: strelka_config
            dbsnp_vcf: dv_dbsnp_vcf
            cosmic_vcf: cosmic_vcf
        out:
            [final_vcf]
