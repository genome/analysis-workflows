#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "running cellranger mkfastq and count"

inputs:
    bcl_directory:
        type: string
    chemistry:
        type: string?
    reference:
        type: string
    sample_name:
        type: string
    simple_sample_csv:
        type: string
    unique_id:
        type: string

steps:
    mkfastq:
        run: ../tools/cellranger_mkfastq.cwl
        in:
            bcl_directory: bcl_directory
            simple_sample_csv: simple_sample_csv
            unique_id: unique_id
        out: [samplesheet_csv, fastq_dir]
    count:
        run: ../tools/cellranger_count.cwl
        in: 
            chemistry: chemistry
            fastq_directory: mkfastq/fastq_dir
            reference: reference
            sample_name: sample_name
            unique_id: unique_id
        out: [out_dir]

outputs:
    counts_out_dir:
        type: Directory
        outputSource: count/out_dir
    fastqs:
        type: Directory
        outputSource: mkfastq/fastq_dir
    samplesheet:
        type: File
        outputSource: mkfastq/samplesheet_csv
