#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "scatter GATK HaplotypeCaller over intervals"
requirements:
    - class: ScatterFeatureRequirement
inputs:
    reference:
        type: string
    cram:
        type: File
        secondaryFiles: [^.crai]
    emit_reference_confidence:
        type: string
    gvcf_gq_bands:
        type: string[]
    intervals:
        type:
            type: array
            items:
                type: array
                items: string

    dbsnp_vcf:
        type: File?
        secondaryFiles: [.tbi]
outputs:
    gvcf:
        type: File[]
        outputSource: haplotype_caller/gvcf
        secondaryFiles: [.tbi]
steps:
    haplotype_caller:
        scatter: [intervals]
        run: gatk_haplotypecaller.cwl
        in:
            reference: reference
            cram: cram
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            dbsnp_vcf: dbsnp_vcf
        out:
            [gvcf]
