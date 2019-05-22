#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Running cellranger count and lineage"

inputs:
    chemistry:
        type: string?
    reference:
        type: Directory
    sample_name:
        type: string
    fastq_directory:
        type: Directory[]
    lineage_min_cell:
        type: int?
        default: 3
    lineage_min_features:
        type: int?
        default: 10

steps:
    count:
        run: ../tools/cellranger_count.cwl
        in:
            chemistry: chemistry
            fastq_directory: fastq_directory
            reference: reference
            sample_name: sample_name
        out: [out_dir]
    lineage:
        run: ../tools/cellmatch_lineage.cwl
        in:
            sample_name: sample_name
            cellranger_out_dir: count/out_dir
            lineage_min_cell: lineage_min_cell
            lineage_min_features: lineage_min_features
        out: [cellmatch_out_dir]

outputs:
    counts_out_dir:
        type: Directory
        outputSource: count/out_dir
    lineage_cellmatch_lineage_out_dir:
        type: Directory
        outputSource: lineage/cellmatch_out_dir
