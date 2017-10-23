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
        type: string
    dbsnp:
        type: File
        secondaryFiles: [.tbi]
    mills:
        type: File
        secondaryFiles: [.tbi]
    known_indels:
        type: File
        secondaryFiles: [.tbi]
outputs:
    final_cram:
        type: File
        outputSource: index_cram/indexed_cram
        secondaryFiles: [.crai, ^.crai]
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
            [tagged_bam]
    merge:
        run: merge.cwl
        in:
            bams: align/tagged_bam
        out:
            [merged_bam]
    position_sort:
        run: position_sort.cwl
        in:
            bam: merge/merged_bam
        out:
            [position_sorted_bam]
    index_bam:
        run: ../detect_variants/index_bam.cwl
        in:
            bam: position_sort/position_sorted_bam
        out:
            [indexed_bam]
    bqsr:
        run: bqsr.cwl
        in:
            reference: reference
            bam: index_bam/indexed_bam
            known_sites: [dbsnp, mills, known_indels]
        out:
            [bqsr_table]
    apply_bqsr:
        run: apply_bqsr.cwl
        in:
            reference: reference
            bam: index_bam/indexed_bam
            bqsr_table: bqsr/bqsr_table
        out:
            [bqsr_bam]
    bam_to_cram:
        run: bam_to_cram.cwl
        in:
            reference: reference
            bam: apply_bqsr/bqsr_bam
        out:
            [cram]
    index_cram:
        run: index_cram.cwl
        in:
            cram: bam_to_cram/cram
        out:
            [indexed_cram]
