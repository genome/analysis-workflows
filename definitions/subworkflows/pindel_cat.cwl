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
    tumor_cram:
        type: File
        secondaryFiles: ["^.crai"]
    normal_cram:
        type: File
        secondaryFiles: ["^.crai"]
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
    tumor_cram_to_bam:
        run: cram_to_bam_and_index.cwl
        in:
            reference: reference
            cram: tumor_cram
        out:
            [bam]
    normal_cram_to_bam:
        run: cram_to_bam_and_index.cwl
        in:
            reference: reference
            cram: normal_cram
        out:
            [bam]
    pindel:
        run: ../tools/pindel.cwl
        in:
            reference: reference
            tumor_bam: tumor_cram_to_bam/bam
            normal_bam: normal_cram_to_bam/bam
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
