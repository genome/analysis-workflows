#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "varscan somatic workflow"
inputs:
    reference:
        type: string
    tumor_cram:
        type: File
        secondaryFiles: [^.crai]
    normal_cram:
        type: File
        secondaryFiles: [^.crai]
    roi_bed:
        type: File?
    strand_filter:
        type: int?
    min_coverage:
        type: int?
    min_var_freq:
        type: float?
    p_value:
        type: float?
    max_normal_freq:
        type: float?
outputs:
    snvs:
        type: File
        outputSource: somatic/snvs
    indels:
        type: File
        outputSource: somatic/indels
    somatic_hc_snvs:
        type: File
        outputSource: process_somatic_snvs/somatic_hc
    somatic_snvs:
        type: File
        outputSource: process_somatic_snvs/somatic
    germline_hc_snvs:
        type: File
        outputSource: process_somatic_snvs/germline_hc
    germline_snvs:
        type: File
        outputSource: process_somatic_snvs/germline
    loh_hc_snvs:
        type: File
        outputSource: process_somatic_snvs/loh_hc
    loh_snvs:
        type: File
        outputSource: process_somatic_snvs/loh
    somatic_hc_indels:
        type: File
        outputSource: process_somatic_indels/somatic_hc
    somatic_indels:
        type: File
        outputSource: process_somatic_indels/somatic
    germline_hc_indels:
        type: File
        outputSource: process_somatic_indels/germline_hc
    germline_indels:
        type: File
        outputSource: process_somatic_indels/germline
    loh_hc_indels:
        type: File
        outputSource: process_somatic_indels/loh_hc
    loh_indels:
        type: File
        outputSource: process_somatic_indels/loh
steps:
    somatic:
        run: ../tools/varscan_somatic.cwl
        in:
            reference: reference
            normal_cram: normal_cram
            tumor_cram: tumor_cram
            roi_bed: roi_bed
            strand_filter: strand_filter
            min_coverage: min_coverage
            min_var_freq: min_var_freq
            p_value: p_value
        out:
            [snvs, indels]
    process_somatic_snvs:
        run: ../tools/varscan_process_somatic.cwl
        in:
            variants: somatic/snvs
            max_normal_freq: max_normal_freq
        out:
            [somatic_hc, somatic, germline_hc, germline, loh_hc, loh]
    process_somatic_indels:
        run: ../tools/varscan_process_somatic.cwl
        in:
            variants: somatic/indels
            max_normal_freq: max_normal_freq
        out:
            [somatic_hc, somatic, germline_hc, germline, loh_hc, loh]
