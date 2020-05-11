#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "given multiple somatic vcfs from the same individual, merges the vcfs, adds readcounts, and creates a table"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: InlineJavascriptRequirement
inputs:
    vcfs:
        type: File[]
        secondaryFiles: [.tbi]
    tumor_bams:
        type: File[]
        secondaryFiles: ${if (self.nameext === ".bam") {return self.basename + ".bai"} else {return self.basename + ".crai"}}
        doc: array of tumor bams or crams
    tumor_sample_names:
        type: string[]
        doc: tumor sample names - assumes same order as the tumor bams
    normal_bam:
        type: File
        secondaryFiles: ${if (self.nameext === ".bam") {return self.basename + ".bai"} else {return self.basename + ".crai"}}
        doc: shared normal sequence file accepts either bam or cram
    normal_sample_name:
        type: string
        doc: shared normal sample name
    reference_fasta:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    min_base_quality:
        type: int?
    min_mapping_quality:
        type: int?
    data_type:
        type:
            - type: enum
              symbols: ["DNA", "RNA"]
        doc: for now, this only accepts either "DNA" or "RNA" and assumes it applies to all samples/bams to avoid having to pass in an array
outputs:
    merged_readcount_vcf:
        type: File
        outputSource: add_readcounts/readcount_vcf
steps:
    merge_vcfs:
        run: ../tools/merge_somatic_vcfs.cwl
        in:
            vcfs: vcfs
            sample_names: tumor_sample_names
            normal_sample_name: normal_sample_name
        out:
            [merged_vcf]

    prepend_normal_bam:
        run: ../tools/prepend_file_to_array.cwl
        in:
            file: normal_bam
            array: tumor_bams
        out:
            [file_array]

    prepend_normal_name:
        run: ../tools/prepend_string_to_array.cwl
        in:
            string: normal_sample_name
            array: tumor_sample_names
        out:
            [string_array]

    add_readcounts:
        run: bam_readcount_multisample.cwl
        in:
            vcf: merge_vcfs/merged_vcf
            bams:
                source: prepend_normal_bam/file_array
            sample_names: 
                source: prepend_normal_name/string_array
            reference_fasta: reference_fasta
            min_base_quality: min_base_quality
            min_mapping_quality: min_mapping_quality
            data_type: data_type
        out:
            [readcount_vcf]
##todo - index the vcf
