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
    bqsr_intervals:
        type: string[]?
    reference:
        type: string
    dbsnp_vcf:
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
    mark_duplicates_metrics_file:
        type: File
        outputSource: mark_duplicates_and_sort/metrics_file
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
    name_sort:
        run: name_sort.cwl
        in:
            bam: merge/merged_bam
        out:
            [name_sorted_bam]
    mark_duplicates_and_sort:
        run: mark_duplicates_and_sort.cwl
        in:
            bam: name_sort/name_sorted_bam
        out:
            [sorted_bam, metrics_file]
    bqsr:
        run: bqsr.cwl
        in:
            reference: reference
            bam: mark_duplicates_and_sort/sorted_bam
            intervals: bqsr_intervals
            known_sites: [dbsnp_vcf, mills, known_indels]
        out:
            [bqsr_table]
    apply_bqsr:
        run: apply_bqsr.cwl
        in:
            reference: reference
            bam: mark_duplicates_and_sort/sorted_bam
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
