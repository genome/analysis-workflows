#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: Workflow
label: "Per-interval pindel"
inputs:
    reference:
        type: File
        secondaryFiles: [".fai"]
    tumor_bam:
        type: File
        secondaryFiles: .bai
    normal_bam:
        type: File
        secondaryFiles: .bai
    interval_list:
        type: File
    insert_size:
        type: int
        default: 400
outputs:
    per_interval_pindel_out:
        type: File
        outputSource: cat/pindel_out
steps:
    sed:
        run: sed.cwl
        in:
           interval_list: interval_list
        out:
           [interval_bed]
    pindel:
        run: pindel.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            insert_size: insert_size
            interval_list: sed/interval_bed
        out:
            [deletions, insertions, tandems, long_insertions, inversions]
    cat:
        run: cat_out.cwl
        in:
            pindel_outs: [pindel/deletions, pindel/insertions, pindel/tandems, pindel/long_insertions, pindel/inversions]
        out:
            [pindel_out]
