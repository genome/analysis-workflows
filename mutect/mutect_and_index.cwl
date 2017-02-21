#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: Workflow
label: "Per-interval mutect"
inputs:
    reference:
        type: File
        secondaryFiles: [".fai", "^.dict"]
    tumor_bam:
        type: File
        secondaryFiles: [^.bai]
    normal_bam:
        type: File?
        secondaryFiles: [^.bai]
    interval_list:
        type: File
    dbsnp_vcf:
        type: File?
        secondaryFiles: [.tbi]
    cosmic_vcf:
        type: File?
        secondaryFiles: [.tbi]
    artifact_detection_mode:
        type: boolean?

outputs:
    vcf:
        type: File
        outputSource: index/indexed_vcf
steps:
    mutect:
        run: mutect.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            interval_list: interval_list
            dbsnp_vcf: dbsnp_vcf
            cosmic_vcf: cosmic_vcf
            artifact_detection_mode: artifact_detection_mode
        out:
            [vcf]
    index:
        run: ../detect_variants/index.cwl
        in:
            vcf: [mutect/vcf]
        out:
            [indexed_vcf]

