#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "varscan somatic workflow"
inputs:
    reference:
        type: File
        secondaryFiles: .fai
    tumor_bam:
        type: File
    normal_bam:
        type: File
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
    mpileup:
        run: samtools_mpileup.cwl
        in:
            reference: reference
            normal_bam: normal_bam
            tumor_bam: tumor_bam
        out:
            [mpileup]
    somatic:
        run: varscan_somatic.cwl
        in:
            mpileup: mpileup/mpileup
        out:
            [snvs, indels]
    process_somatic_snvs:
        run: varscan_processsomatic.cwl
        in:
            variants: somatic/snvs
        out:
            [somatic_hc, somatic, germline_hc, germline, loh_hc, loh]
    process_somatic_indels:
        run: varscan_processsomatic.cwl
        in:
            variants: somatic/indels
        out:
            [somatic_hc, somatic, germline_hc, germline, loh_hc, loh]
