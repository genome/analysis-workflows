#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "HS Metrics workflow"
requirements:
    - class: ScatterFeatureRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/labelled_file.yml
    - class: StepInputExpressionRequirement
    - class: SubworkflowFeatureRequirement
inputs:
    cram:
        type: File
        secondaryFiles: [^.crai]
    minimum_mapping_quality:
        type: int?
    minimum_base_quality:
        type: int?
    per_base_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    per_target_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
    reference:
        type: string
    summary_intervals:
        type: ../types/labelled_file.yml#labelled_file[]
outputs:
    per_base_coverage_metrics:
        type: File[]
        outputSource: collect_per_base_hs_metrics/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File[]
        outputSource: collect_per_base_hs_metrics/hs_metrics
    per_target_coverage_metrics:
        type: File[]
        outputSource: collect_per_target_hs_metrics/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File[]
        outputSource: collect_per_target_hs_metrics/hs_metrics
    summary_hs_metrics:
        type: File[]
        outputSource: collect_summary_hs_metrics/hs_metrics
steps:
    collect_summary_hs_metrics:
        run: ../tools/collect_hs_metrics.cwl
        scatter: [bait_intervals, target_intervals, output_prefix]
        scatterMethod: dotproduct
        in:
            cram: cram
            reference: reference
            metric_accumulation_level:
                valueFrom: "ALL_READS"
            bait_intervals:
                source: summary_intervals
                valueFrom: $(self.file)
            target_intervals:
                source: summary_intervals
                valueFrom: $(self.file)
            per_target_coverage:
                default: false
            per_base_coverage:
                default: false
            output_prefix:
                source: summary_intervals
                valueFrom: "summary-$(self.label)"
            minimum_mapping_quality: minimum_mapping_quality
            minimum_base_quality: minimum_base_quality
        out:
            [hs_metrics]
    collect_per_base_hs_metrics:
        run: ../tools/collect_hs_metrics.cwl
        scatter: [bait_intervals, target_intervals, output_prefix]
        scatterMethod: dotproduct
        in:
            cram: cram
            reference: reference
            metric_accumulation_level:
                valueFrom: "ALL_READS"
            bait_intervals:
                source: per_base_intervals
                valueFrom: $(self.file)
            target_intervals:
                source: per_base_intervals
                valueFrom: $(self.file)
            per_target_coverage:
                default: false
            per_base_coverage:
                default: true
            output_prefix:
                source: per_base_intervals
                valueFrom: "base-$(self.label)"
            minimum_mapping_quality: minimum_mapping_quality
            minimum_base_quality: minimum_base_quality
        out:
            [hs_metrics, per_base_coverage_metrics]
    collect_per_target_hs_metrics:
        run: ../tools/collect_hs_metrics.cwl
        scatter: [bait_intervals, target_intervals, output_prefix]
        scatterMethod: dotproduct
        in:
            cram: cram
            reference: reference
            metric_accumulation_level:
                valueFrom: "ALL_READS"
            bait_intervals:
                source: per_target_intervals
                valueFrom: $(self.file)
            target_intervals:
                source: per_target_intervals
                valueFrom: $(self.file)
            per_target_coverage:
                default: true
            per_base_coverage:
                default: false
            output_prefix:
                source: per_target_intervals
                valueFrom: "target-$(self.label)"
            minimum_mapping_quality: minimum_mapping_quality
            minimum_base_quality: minimum_base_quality
        out:
            [hs_metrics, per_target_coverage_metrics]
