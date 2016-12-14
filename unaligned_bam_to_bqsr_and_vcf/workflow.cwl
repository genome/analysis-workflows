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
        secondaryFiles: [.fai, .bwt, .sa, .ann, .amb, .pac, ^.dict]
    tmpdir:
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
    final_bam:
        type: File
        outputSource: apply_bqsr/bqsr_bam
        secondaryFiles: [.bai]
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
            bams: [align/tagged_bam]
        out:
            [merged_bam]
    name_sort:
        run: name_sort.cwl
        in:
            bam: merge/merged_bam
            tmpdir: tmpdir
        out:
            [name_sorted_bam]
    mark_duplicates:
        run: mark_duplicates.cwl
        in:
            bam: name_sort/name_sorted_bam
        out:
            [duplicate_marked_bam]
    sort:
        run: sort.cwl
        in:
            bam: mark_duplicates/duplicate_marked_bam
            tmpdir: tmpdir
        out:
            [sorted_bam]
    bqsr:
        run: bqsr.cwl
        in:
            reference: reference
            bam: sort/sorted_bam
            known_sites: [dbsnp, mills, known_indels]
        out:
            [bqsr_table]
    apply_bqsr:
        run: apply_bqsr.cwl
        in:
            reference: reference
            bam: sort/sorted_bam
            bqsr_table: bqsr/bqsr_table
        out:
            [bqsr_bam]
