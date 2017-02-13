#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: Workflow
label: "Per-chromosome pindel"
inputs:
    reference:
        type: File
        secondaryFiles: [".fai"]
    tumor_bam:
        type: File
        secondaryFiles: ["^.bai"]
    normal_bam:
        type: File
        secondaryFiles: ["^.bai"]
    chromosome:
        type: string
    insert_size:
        type: int
        default: 400
outputs:
    per_chromosome_pindel_out:
        type: File
        outputSource: cat/pindel_out
steps:
    pindel:
        run: pindel.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            insert_size: insert_size
            chromosome: chromosome
        out:
            [deletions, insertions, tandems, long_insertions, inversions]
    cat:
        run: cat_out.cwl
        in:
            pindel_outs: [pindel/deletions, pindel/insertions, pindel/tandems, pindel/long_insertions, pindel/inversions]
        out:
            [pindel_out]
