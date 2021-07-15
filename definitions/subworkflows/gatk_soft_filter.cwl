#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "apply soft filtering to a gatk called vcf using hard filter paramaters"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: MultipleInputFeatureRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    vcf:
        type: File
        secondaryFiles: [.tbi]
outputs:
    filtered_vcf:
        type: File
        secondaryFiles: [.tbi]
        outputSource: index_merged/indexed_vcf
steps:
    split_snps:
        run: ../tools/select_variants.cwl
        in:
            reference: reference
            vcf: vcf
            output_vcf_basename:
                default: "SNPS"
            select_type:
                default: "SNP"
        out:
            [filtered_vcf]
    split_indels:
        run: ../tools/select_variants.cwl
        in:
            reference: reference
            vcf: vcf
            output_vcf_basename:
                default: "INDELS"
            select_type:
                default: "INDEL"
        out:
            [filtered_vcf]
    filter_snps:
        run: ../tools/variant_filtration.cwl
        in:
            reference: reference
            vcf: split_snps/filtered_vcf
            filters:
                default: ["QD<2.0;QD2", "QUAL<30.0;QUAL30", "SOR>3.0;SOR3", "FS>60.0;FS60", "MQ<40.0;MQ40", "MQRankSum<-12.5;MQRankSum-12.5", "ReadPosRankSum<-8.0;ReadPosRankSum-8"]
            output_vcf_basename:
                default: "SNPS.filtered"
        out:
            [filtered_vcf]
    filter_indels:
        run: ../tools/variant_filtration.cwl
        in:
            reference: reference
            vcf: split_indels/filtered_vcf
            filters:
                default: ["QD<2.0;QD2", "QUAL<30.0;QUAL30", "FS>200.0;FS200", "ReadPosRankSum<-20.0;ReadPosRankSum-20"]
            output_vcf_basename:
                default: "INDELS.filtered"
        out:
            [filtered_vcf]
    merge:
        run: ../tools/merge_vcf.cwl
        in:
            merged_vcf_basename:
                default: "soft_filtered"
            vcfs:
                source: [filter_snps/filtered_vcf, filter_indels/filtered_vcf]
                linkMerge: merge_flattened
        out:
            [merged_vcf]
    index_merged:
        run: ../tools/index_vcf.cwl
        in:
            vcf: merge/merged_vcf
        out:
            [indexed_vcf]
