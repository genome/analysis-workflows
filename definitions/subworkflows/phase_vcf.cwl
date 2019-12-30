#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "phase VCF"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement
inputs:
    somatic_vcf:
        type: File
        secondaryFiles: [.tbi]
    germline_vcf:
        type: File
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    reference_dict:
        type: File
    bam:
        type: File
        secondaryFiles: ${if (self.nameext === ".bam") {return self.basename + ".bai"} else {return self.basename + ".crai"}}
    normal_sample_name:
        type: string
    tumor_sample_name:
        type: string
outputs:
    phased_vcf:
        type: File
        outputSource: bgzip_and_index_phased_vcf/indexed_vcf
        secondaryFiles: [.tbi]
steps:
    rename_germline_vcf:
        run: ../tools/replace_vcf_sample_name.cwl
        in:
            input_vcf: germline_vcf
            sample_to_replace: normal_sample_name
            new_sample_name: tumor_sample_name
        out: [renamed_vcf]
    index_renamed_germline:
        run: ../tools/index_vcf.cwl
        in:
            vcf: rename_germline_vcf/renamed_vcf
        out:
            [indexed_vcf]

    select_somatic_tumor_sample:
        run: ../tools/select_variants.cwl
        in:
            reference: reference
            vcf: somatic_vcf
            output_vcf_basename:
                default: 'somatic_tumor_only'
            samples_to_include:
                source: tumor_sample_name
                valueFrom: ${ return [self]; }
        out:
            [filtered_vcf]
    index_filtered_somatic:
        run: ../tools/index_vcf.cwl
        in:
            vcf: select_somatic_tumor_sample/filtered_vcf
        out: [indexed_vcf]
    combine_variants:
        run: ../tools/pvacseq_combine_variants.cwl
        in:
            reference: reference
            germline_vcf: index_renamed_germline/indexed_vcf
            somatic_vcf: index_filtered_somatic/indexed_vcf
        out:
            [combined_vcf]
    sort:
        run: ../tools/sort_vcf.cwl
        in:
            vcf: combine_variants/combined_vcf
            reference_dict: reference_dict
        out:
            [sorted_vcf]
    bgzip_and_index:
        run: bgzip_and_index.cwl
        in:
            vcf: sort/sorted_vcf
        out:
            [indexed_vcf]
    phase_vcf:
        run: ../tools/read_backed_phasing.cwl
        in:
            reference: reference
            bam: bam
            vcf: bgzip_and_index/indexed_vcf
        out:
            [phased_vcf]
    bgzip_and_index_phased_vcf:
        run: bgzip_and_index.cwl
        in:
            vcf: phase_vcf/phased_vcf
        out:
            [indexed_vcf]
