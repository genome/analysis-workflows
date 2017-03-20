#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: Workflow
label: "Per-chromosome pindel"
requirements:
    - class: MultipleInputFeatureRequirement
inputs:
    reference:
        type: File
        secondaryFiles: [".fai"]
    tumor_bam:
        type: File
    normal_bam:
        type: File
    tumor_bam_index:
        type: File
    normal_bam_index:
        type: File
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
            tumor_bam_index: tumor_bam_index
            normal_bam_index: normal_bam_index
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
