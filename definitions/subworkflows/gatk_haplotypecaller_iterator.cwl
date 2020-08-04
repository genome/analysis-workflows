#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "scatter GATK HaplotypeCaller over intervals"
requirements:
    - class: InlineJavascriptRequirement
    - class: ScatterFeatureRequirement
    - class: StepInputExpressionRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    bam:
        type: File
        secondaryFiles: [^.bai]
    emit_reference_confidence:
        type:
            type: enum
            symbols: ['NONE', 'BP_RESOLUTION', 'GVCF']
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
    contamination_fraction:
        type: string?
    max_alternate_alleles:
        type: int?
    ploidy:
        type: int?
    read_filter:
        type: string?
    output_prefix:
        type: string?
        doc: 'an optional string to prepend to the output filename'
outputs:
    gvcf:
        type: File[]
        outputSource: haplotype_caller/gvcf
        secondaryFiles: [.tbi]
steps:
    haplotype_caller:
        scatter: [intervals]
        run: ../tools/gatk_haplotype_caller.cwl
        in:
            reference: reference
            bam: bam
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            dbsnp_vcf: dbsnp_vcf
            contamination_fraction: contamination_fraction
            max_alternate_alleles: max_alternate_alleles
            ploidy: ploidy
            read_filter: read_filter
            output_file_name:
                source: output_prefix
                valueFrom: '${
                    var prefix = self !== null ? self : "";
                    if (inputs.intervals.length == 1 && inputs.intervals[0].match(/^[0-9A-Za-z]+$/)) {
                        return prefix + inputs.intervals[0] + ".g.vcf.gz";
                    } else {
                        return prefix + "output.g.vcf.gz";
                    }
                }'
        out:
            [gvcf]
