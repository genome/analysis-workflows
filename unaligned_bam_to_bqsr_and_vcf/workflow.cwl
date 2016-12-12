#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Unaligned BAM to BQSR and VCF"
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement

inputs:
    bams:
        type: File[]
    readgroups:
        type: string[]
    reference:
        type: File
        secondaryFiles: [.fai, .bwt, .sa, .ann, .amb, .pac]
outputs:
    merged_bam:
        type: File
        outputSource: merge/merged_bam
steps:
    align:
        scatter: [bam, readgroup]
        scatterMethod: dotproduct
        run: align.cwl
        in:
            bam: bams
            readgroup: readgroups
            reference: reference
        out:
            [tagged_sam]
    merge:
        run: merge.cwl
        in:
            bams: [align/tagged_sam]
        out:
            [merged_bam]
