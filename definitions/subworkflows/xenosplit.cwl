#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Xenosplit workflow"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    graftbam:
        type: File
        secondaryFiles: [.bai]
    hostbam:
        type: File
        secondaryFiles: [.bai]
outputs:
    xenosplitbam:
        type: File
        outputSource: xenosplit/graftOut
    xenosplitscore:
        type: File
        outputSource: xenosplit/goodnessOfMapping
steps:
    xenosplit_bam_conversion:
        run: ../tools/xenosplit_bam_conversion.cwl
        in:
            graftbam: graftbam
            hostbam: hostbam
        out:
            [graftbam_accepted, hostbam_accepted]
    xenosplit:
        run: ../tools/xenosplit.cwl
        in:
            graftbam: xenosplit_bam_conversion/graftbam_accepted
            hostbam: xenosplit_bam_conversion/hostbam_accepted
        out:
            [graftOut, goodnessOfMapping]
