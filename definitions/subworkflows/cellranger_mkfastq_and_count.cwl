#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "running cellranger mkfastq and count"

inputs:
    bcl_directory:
        type: Directory
    chemistry:
        type: string?
    reference:
        type: Directory
    sample_name:
        type: string
    simple_sample_csv:
        type: File

steps:
    mkfastq:
        run: ../tools/cellranger_mkfastq.cwl
        in:
            bcl_directory: bcl_directory
            simple_sample_csv: simple_sample_csv
        out: [samplesheet_csv, fastq_dir]
    count:
        run: ../tools/cellranger_count.cwl
        in: 
            chemistry: chemistry
            fastq_directory: mkfastq/fastq_dir
            reference: reference
            sample_name: sample_name
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
