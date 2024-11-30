#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "filter jointly called vcfs from read based callers"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement
    - class: ScatterFeatureRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    sample_names:
        type: string[]
    bams:
        type: File[]
        secondaryFiles: [^.bai]
    filter_del_depth:
        type: double?
    filter_dup_depth:
        type: double?
    filter_paired_count:
        type: int?
    filter_split_count:
        type: int?
    filter_alt_abundance_percentage:
        type: double?
    sv_vcf:
        type: File
        secondaryFiles: [.tbi]
    vcf_source:
        type:
          - type: enum
            symbols: ["manta", "smoove"]
outputs:
    vcfs:
        type: File[]
        outputSource: final_index/indexed_vcf
        secondaryFiles: [.tbi]
steps:
    read_support_filter:
        run: ../tools/filter_sv_vcf_read_support.cwl
        in:
            abundance_percentage: filter_alt_abundance_percentage
            input_vcf: sv_vcf
            paired_count: filter_paired_count
            split_count: filter_split_count
            vcf_source: vcf_source
        out:
            [filtered_sv_vcf]
    bgzip_index:
        run: bgzip_and_index.cwl
        in:
            vcf: read_support_filter/filtered_sv_vcf
        out:
            [indexed_vcf]
    split_vcf:
        scatter: [sample_name]
        run: ../tools/bcftools_view.cwl
        in:
            sample_name: sample_names
            in_vcf: bgzip_index/indexed_vcf
        out:
            [vcf]
    duphold:
        scatter: [bam, sv_vcf]
        scatterMethod: dotproduct
        run: ../tools/duphold.cwl
        in:
            bam: bams
            reference: reference
            sv_vcf: split_vcf/vcf
        out:
            [annotated_sv_vcf]
    bgzip_index_duphold:
        scatter: [vcf]
        scatterMethod: dotproduct
        run: bgzip_and_index.cwl
        in:
            vcf: duphold/annotated_sv_vcf
        out:
            [indexed_vcf]
    merge_vcfs:
        run: ../tools/bcftools_merge.cwl
        in:
            vcfs: bgzip_index_duphold/indexed_vcf
        out:
            [merged_vcf]
    depth_filter:
        run: ../tools/filter_sv_vcf_depth.cwl
        in:
            input_vcf: merge_vcfs/merged_vcf
            deletion_depth: filter_del_depth
            duplication_depth: filter_dup_depth
            vcf_source:
                default: "duphold"
        out:
            [filtered_sv_vcf]
    final_split_vcf:
        scatter: [sample_name, output_vcf_name]
        scatterMethod: dotproduct
        run: ../tools/bcftools_view.cwl
        in:
            sample_name: sample_names
            in_vcf: depth_filter/filtered_sv_vcf
            vcf_source: vcf_source
            output_vcf_name:
                source: [sample_names]
                valueFrom: |
                    ${
                      var sample = self;
                      var caller = inputs.vcf_source;
                      var result = sample + "-" + caller + ".vcf.gz";
                      return result;
                    }
        out:
            [vcf]
    rename:
        scatter: [input_vcf, sample_to_replace, new_sample_name, output_name]
        scatterMethod: dotproduct
        run: ../tools/replace_vcf_sample_name.cwl
        in:
            input_vcf: final_split_vcf/vcf
            sample_to_replace: sample_names
            vcf_source: vcf_source
            new_sample_name:
                source: [sample_names]
                valueFrom: |
                    ${
                      var sample = self;
                      var caller = inputs.vcf_source;
                      var result = sample + "-" + caller;
                      return result;
                    }
            output_name:
                source: [sample_names]
                valueFrom: |
                    ${
                      var sample = self;
                      var caller = inputs.vcf_source;
                      var result = sample + "-" + caller + ".vcf.gz";
                      return result;
                    }
        out:
            [renamed_vcf]
    final_index:
        scatter: [vcf]
        scatterMethod: dotproduct
        run: ../tools/index_vcf.cwl
        in:
            vcf: rename/renamed_vcf
        out:
            [indexed_vcf]
