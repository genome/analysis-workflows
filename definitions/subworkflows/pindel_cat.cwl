#!/usr/bin/env cwl-runner

cwlVersion: v1.0

class: Workflow
label: "Per-region pindel"
requirements:
    - class: MultipleInputFeatureRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    tumor_bam:
        type: File
        secondaryFiles: ["^.bai"]
    normal_bam:
        type: File
        secondaryFiles: ["^.bai"]
    region_file:
        type: File
    insert_size:
        type: int
        default: 400
    tumor_sample_name:
        type: string
    normal_sample_name:
        type: string
outputs:
    per_region_pindel_out:
        type: File
        outputSource: cat/pindel_out
steps:
    pindel:
        run: ../tools/pindel.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            insert_size: insert_size
            region_file: region_file
            tumor_sample_name: tumor_sample_name
            normal_sample_name: normal_sample_name
        out:
            [deletions, insertions, tandems, long_insertions, inversions]
    cat:
        run: ../tools/cat_out.cwl
        in:
            pindel_outs: [pindel/deletions, pindel/insertions, pindel/tandems, pindel/long_insertions, pindel/inversions]
        out:
            [pindel_out]
