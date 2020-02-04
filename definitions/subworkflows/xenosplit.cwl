#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Xenosplit workflow"
requirements:
    - class: SubworkflowFeatureRequirement

inputs:
    graft_reference:
        type: Directory
    host_reference:
        type: Directory
    fastq:
        type: File[]
    fastq2:
        type: File[]
    graft_outfile_name_prefix:
        type: string?
        default: "Graft_"
    host_outfile_name_prefix:
        type: string?
        default: "Host_"

outputs:
    xenosplitbam:
        type: File
        outputSource: xenosplit/graftOut
    xenosplitscore:
        type: File
        outputSource: xenosplit/goodnessOfMapping

steps:
    graft_star_alignment:
        run: ../tools/xenosplit_align.cwl
        in:
            fastq: fastq
            fastq2: fastq2
            star_genome_dir: graft_reference
            outfile_name_prefix: graft_outfile_name_prefix
        out:
            [aligned_bam]
    host_star_alignment:
        run: ../tools/xenosplit_align.cwl
        in:
            fastq: fastq
            fastq2: fastq2
            star_genome_dir: host_reference
            outfile_name_prefix: host_outfile_name_prefix
        out:
            [aligned_bam]
    xenosplit_bam_conversion:
        run: ../tools/xenosplit_bam_conversion.cwl
        in:
            graftbam: graft_star_alignment/aligned_bam
            hostbam: host_star_alignment/aligned_bam
        out:
            [graftbam_accepted, hostbam_accepted]
    xenosplit:
        run: ../tools/xenosplit.cwl
        in:
            graftbam: xenosplit_bam_conversion/graftbam_accepted
            hostbam: xenosplit_bam_conversion/hostbam_accepted
        out:
            [graftOut, goodnessOfMapping]
