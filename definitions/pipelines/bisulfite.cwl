#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Bisulfite alignment and QC"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
          - $import: ../types/trimming_options.yml
inputs:
    reference_index:
        type: string
    reference_sizes:
        type: File
    sequence:
        type: ../types/sequence_data.yml#sequence_data[]
        doc: |
          sequence represents the sequencing data as either FASTQs or BAMs with accompanying
          readgroup information. Note that in the @RG field ID and SM are required.
    sample_name:
        type: string
    trimming_options:
        type:
            - ../types/trimming_options.yml#trimming_options
            - "null"
    QCannotation:
        type: File
    assay_non_cpg_sites:
        type:
            type: enum
            symbols: ["true", "false"]
        default: "false"
        doc: "Variable to determine if user wants to obtain obtain bed/bigwig files for non-CpG cytosines. Value - true or false"
outputs:
    cram:
        type: File
        outputSource: index_cram/indexed_cram
        secondaryFiles: [.crai, ^.crai]
    vcf:
        type: File
        outputSource: pileup/vcf
    cpgs:
        type: File[]
        outputSource: vcf2bed/methylation_bed
    cpg_bigwig:
        type: File[]
        outputSource: bedgraph_to_bigwig/methylation_bigwig
    gathered_directory:
        type: Directory
        outputSource: bisulfite_qc/QC_directory
steps:
    bisulfite_alignment:
        run: ../subworkflows/sequence_to_bisulfite_alignment.cwl
        in:
            sequence: sequence
            trimming_options: trimming_options
            reference_index: reference_index
            sample_name: sample_name
        out:
            [aligned_bam]
    pileup:
        run: ../tools/biscuit_pileup.cwl
        in:
            bam: bisulfite_alignment/aligned_bam
            reference: reference_index
        out:
            [vcf]
    bisulfite_qc:
        run: ../subworkflows/bisulfite_qc.cwl
        in:
            vcf: pileup/vcf
            bam: bisulfite_alignment/aligned_bam
            reference: reference_index
            QCannotation: QCannotation
        out:
            [QC_directory]
    vcf2bed:
        run: ../tools/bisulfite_vcf2bed.cwl
        in:
            vcf: pileup/vcf
            reference: reference_index
            assay_non_cpg_sites: assay_non_cpg_sites
        out:
            [methylation_bed,methylation_bedgraph]
    bedgraph_to_bigwig:
        run: ../tools/bedgraph_to_bigwig.cwl
        in:
            methylation_bedgraph: vcf2bed/methylation_bedgraph
            reference_sizes: reference_sizes
        out:
            [methylation_bigwig]
    bam_to_cram:
        run: ../tools/bam_to_cram.cwl
        in:
            reference: reference_index
            bam: bisulfite_alignment/aligned_bam
        out:
            [cram]
    index_cram:
        run: ../tools/index_cram.cwl
        in:
            cram: bam_to_cram/cram
        out:
            [indexed_cram]
