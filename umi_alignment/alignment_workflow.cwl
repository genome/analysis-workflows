#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "umi per-lane alignment subworkflow"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    bam:
        type: File
    read_structure:
        type: string[]
    reference: string
outputs:
    aligned_bam:
        type: File
        secondaryFiles: [^.bai]
        outputSource: align/aligned_bam
    adapter_metrics:
        type: File
        outputSource: mark_illumina_adapters/metrics
steps:
    extract_umis:
        run: ../definitions/tools/extract_umis.cwl
        in:
            bam: bam
            read_structure: read_structure
        out:
            [umi_extracted_bam]
    mark_illumina_adapters:
        run: mark_illumina_adapters.cwl
        in:
            bam: extract_umis/umi_extracted_bam
        out:
            [marked_bam, metrics]
    align:
        run: align.cwl
        in:
            bam: mark_illumina_adapters/marked_bam
            reference: reference
        out:
            [aligned_bam]
