#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Varscan Workflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: MultipleInputFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    tumor_bam:
        type: File
        secondaryFiles: [^.bai, .bai]
    normal_bam:
        type: File
        secondaryFiles: [^.bai]
    interval_list:
        type: File
    strand_filter:
        type: int?
        default: 0
    min_coverage:
        type: int?
        default: 8
    min_var_freq:
        type: float?
        default: 0.1
    p_value:
        type: float?
        default: 0.99
    max_normal_freq:
        type: float?
    normal_sample_name:
        type: string
    tumor_sample_name:
        type: string
    scatter_count:
        type: int
        default: 50
outputs:
    unfiltered_vcf:
        type: File
        outputSource: filter/unfiltered_vcf
        secondaryFiles: [.tbi]
    filtered_vcf:
        type: File
        outputSource: filter/filtered_vcf
        secondaryFiles: [.tbi]
steps:
    split_interval_list:
        run: ../tools/split_interval_list.cwl
        in:
            interval_list: interval_list
            scatter_count: scatter_count
        out: [split_interval_lists]
    intervals_to_bed:
        scatter: interval_list
        run: ../tools/intervals_to_bed.cwl
        in:
            interval_list: split_interval_list/split_interval_lists
        out:
            [interval_bed]
    varscan:
        scatter: roi_bed
        run: varscan.cwl
        in:
            reference: reference
            tumor_bam: tumor_bam
            normal_bam: normal_bam
            roi_bed: intervals_to_bed/interval_bed
            strand_filter: strand_filter
            min_coverage: min_coverage
            min_var_freq: min_var_freq
            p_value: p_value
            max_normal_freq: max_normal_freq
        out:
            [somatic_snvs, somatic_indels, somatic_hc_snvs, somatic_hc_indels]
    merge_scattered_somatic_snvs:
        run: ../tools/picard_merge_vcfs.cwl
        in:
            vcfs: varscan/somatic_snvs
            sequence_dictionary:
                source: reference
                valueFrom: '$( (self.secondaryFiles !== undefined) ? self.secondaryFiles.find(function(x) { return x.nameext === ".dict" }) : self.replace(/.fa$/,".dict") )'
            merged_vcf_basename:
                valueFrom: 'somatic_snvs'
        out:
            [merged_vcf]
    merge_scattered_somatic_indels:
        run: ../tools/picard_merge_vcfs.cwl
        in:
            vcfs: varscan/somatic_indels
            sequence_dictionary:
                source: reference
                valueFrom: '$((self.secondaryFiles !== undefined) ? self.secondaryFiles.find(function(x) { return x.nameext === ".dict" }) : self.replace(/.fa$/,".dict") )'
            merged_vcf_basename:
                valueFrom: 'somatic_indels'
        out:
            [merged_vcf]
    merge_scattered_somatic_hc_snvs:
        run: ../tools/picard_merge_vcfs.cwl
        in:
            vcfs: varscan/somatic_hc_snvs
            sequence_dictionary:
                source: reference
                valueFrom: '$((self.secondaryFiles !== undefined) ? self.secondaryFiles.find(function(x) { return x.nameext === ".dict" }) : self.replace(/.fa$/,".dict") )'
            merged_vcf_basename:
                valueFrom: 'somatic_hc_snvs'
        out:
            [merged_vcf]
    merge_scattered_somatic_hc_indels:
        run: ../tools/picard_merge_vcfs.cwl
        in:
            vcfs: varscan/somatic_hc_indels
            sequence_dictionary:
                source: reference
                valueFrom: '$((self.secondaryFiles !== undefined) ? self.secondaryFiles.find(function(x) { return x.nameext === ".dict" }) : self.replace(/.fa$/,".dict") )'
            merged_vcf_basename:
                valueFrom: 'somatic_hc_indels'
        out:
            [merged_vcf]
    merge_snvs:
        run: ../tools/set_filter_status.cwl
        in:
            vcf: merge_scattered_somatic_snvs/merged_vcf
            filtered_vcf: merge_scattered_somatic_hc_snvs/merged_vcf
            reference: reference
        out:
            [merged_vcf]
    index_merged_snvs:
        run: ../tools/index_vcf.cwl
        in:
            vcf: merge_snvs/merged_vcf
        out:
            [indexed_vcf]
    merge_indels:
        run: ../tools/set_filter_status.cwl
        in:
            vcf: merge_scattered_somatic_indels/merged_vcf
            filtered_vcf: merge_scattered_somatic_hc_indels/merged_vcf
            reference: reference
        out:
            [merged_vcf]
    index_merged_indels:
        run: ../tools/index_vcf.cwl
        in:
            vcf: merge_indels/merged_vcf
        out:
            [indexed_vcf]
    merge:
        run: ../tools/merge_vcf.cwl
        in:
            vcfs: [index_merged_snvs/indexed_vcf, index_merged_indels/indexed_vcf]
        out:
            [merged_vcf]
    rename_tumor_sample:
        run: ../tools/replace_vcf_sample_name.cwl
        in: 
            input_vcf: merge/merged_vcf
            sample_to_replace:
                default: 'TUMOR'
            new_sample_name: tumor_sample_name
        out:
            [renamed_vcf]
    rename_normal_sample:
        run: ../tools/replace_vcf_sample_name.cwl
        in: 
            input_vcf: rename_tumor_sample/renamed_vcf
            sample_to_replace:
                default: 'NORMAL'
            new_sample_name: normal_sample_name
        out:
            [renamed_vcf]
    index:
        run: ../tools/index_vcf.cwl
        in:
            vcf: rename_normal_sample/renamed_vcf
        out:
            [indexed_vcf]
    filter:
        run: fp_filter.cwl
        in:
            reference: reference
            bam: tumor_bam
            vcf: index/indexed_vcf
            min_var_freq: min_var_freq
            variant_caller: 
                valueFrom: "varscan"
            sample_name: tumor_sample_name
        out:
            [unfiltered_vcf, filtered_vcf]
