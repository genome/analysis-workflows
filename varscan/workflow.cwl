#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Varscan Workflow"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement
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
        outputSource: index_snvs/indexed_vcf
        secondaryFiles: .tbi
    indels:
        type: File
        outputSource: index_indels/indexed_vcf
        secondaryFiles: .tbi
    merged_vcf:
        type: File
        outputSource: index/indexed_vcf
        secondaryFiles: .tbi
steps:
    varscan:
        run: varscan.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
        out:
            [somatic_snvs, somatic_indels, somatic_hc_snvs, somatic_hc_indels]
    bgzip_and_index_snvs:
        run: bgzip_and_index.cwl
        in:
            vcf: varscan/somatic_snvs
        out:
            [indexed_vcf]
    bgzip_and_index_hc_snvs:
        run: bgzip_and_index.cwl
        in:
            vcf: varscan/somatic_hc_snvs
        out:
            [indexed_vcf]
    bgzip_and_index_indels:
        run: bgzip_and_index.cwl
        in:
            vcf: varscan/somatic_indels
        out:
            [indexed_vcf]
    bgzip_and_index_hc_indels:
        run: bgzip_and_index.cwl
        in:
            vcf: varscan/somatic_hc_indels
        out:
            [indexed_vcf]
    merge_snvs:
        run: set_filter_status.cwl
        in:
            vcf: bgzip_and_index_snvs/indexed_vcf
            filtered_vcf: bgzip_and_index_hc_snvs/indexed_vcf
            reference: reference
        out:
            [merged_vcf]
    index_snvs:
        run: ../detect_variants/index.cwl
        in:
            vcf: merge_snvs/merged_vcf
        out:
            [indexed_vcf]
    merge_indels:
        run: set_filter_status.cwl
        in:
            vcf: bgzip_and_index_indels/indexed_vcf
            filtered_vcf: bgzip_and_index_hc_indels/indexed_vcf
            reference: reference
        out:
            [merged_vcf]
    index_indels:
        run: ../detect_variants/index.cwl
        in:
            vcf: merge_indels/merged_vcf
        out:
            [indexed_vcf]
    merge:
        run: ../detect_variants/merge.cwl
        in:
            vcfs: [index_snvs/indexed_vcf, index_indels/indexed_vcf]
        out:
            [merged_vcf]
    index:
        run: ../detect_variants/index.cwl
        in:
            vcf: merge/merged_vcf
        out:
            [indexed_vcf]
