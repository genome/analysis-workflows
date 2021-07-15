#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Filter multiple sv vcfs from depth callers(cnvkit/cnvnator), returns single sample vcfs with the sample name as $SAMPLE-$CALLER"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement
    - class: ScatterFeatureRequirement
inputs:
    bams:
        type: File[]
        secondaryFiles: [^.bai]
    sample_names:
        type: string[]
    filter_del_depth:
        type: double?
    filter_dup_depth:
        type: double?
    min_sv_size:
        type: int?
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    sv_vcfs:
        type: File[]
    vcf_source:
        type:
          - type: enum
            symbols: ["cnvkit", "cnvnator"]
    merge_distance:
        type: int?
outputs:
    vcfs:
        type: File[]
        outputSource: bgzip_and_index/indexed_vcf
        secondaryFiles: [.tbi]
steps:
    merge_calls:
        scatter: [input_vcf]
        run: ../tools/custom_merge_sv_records.cwl
        in:
            input_vcf: sv_vcfs
            distance: merge_distance
        out:
            [vcf]
    size_filter:
        scatter: [input_vcf]
        run: ../tools/filter_sv_vcf_size.cwl
        in:
            input_vcf: merge_calls/vcf
            size_method:
                default: "min_len"
            sv_size: min_sv_size
        out:
            [filtered_sv_vcf]
    duphold:
        scatter: [bam, sv_vcf]
        scatterMethod: dotproduct
        run: ../tools/duphold.cwl
        in:
            bam: bams
            reference: reference
            sv_vcf: size_filter/filtered_sv_vcf
        out:
            [annotated_sv_vcf]
    depth_filter:
        scatter: [input_vcf, output_vcf_name]
        scatterMethod: dotproduct
        run: ../tools/filter_sv_vcf_depth.cwl
        in:
            input_vcf: duphold/annotated_sv_vcf
            deletion_depth: filter_del_depth
            duplication_depth: filter_dup_depth
            output_vcf_name:
                source: [sample_names]
                valueFrom: |
                    ${
                      var sample = self;
                      var caller = inputs.vcf_source;
                      var vcf_name = sample + "-" + caller  + ".vcf";
                      return vcf_name;
                    }
            vcf_source:
                default: "duphold"
        out:
            [filtered_sv_vcf]
    rename:
        scatter: [input_vcf, new_sample_name, sample_to_replace, output_name]
        scatterMethod: dotproduct
        run: ../tools/replace_vcf_sample_name.cwl
        in:
            input_vcf: depth_filter/filtered_sv_vcf
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
    bgzip_and_index:
        scatter: [vcf]
        run: bgzip_and_index.cwl
        in:
            vcf: rename/renamed_vcf
        out: [indexed_vcf]
