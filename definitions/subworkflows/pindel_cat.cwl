#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: Workflow
label: "Per-chromosome pindel"
requirements:
    - class: MultipleInputFeatureRequirement
inputs:
    reference:
        type: string
    tumor_cram:
        type: File
        secondaryFiles: ["^.crai"]
    normal_cram:
        type: File
        secondaryFiles: ["^.crai"]
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
        run: ../tools/pindel.cwl
        in:
            reference: reference
            tumor_cram: tumor_cram
            normal_cram: normal_cram
            insert_size: insert_size
            chromosome: chromosome
        out:
            [deletions, insertions, tandems, long_insertions, inversions]
    cat:
        run: ../tools/cat_out.cwl
        in:
            pindel_outs: [pindel/deletions, pindel/insertions, pindel/tandems, pindel/long_insertions, pindel/inversions]
        out:
            [pindel_out]
