#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Downsample and HaplotypeCaller"
requirements:
    - class: ScatterFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: SubworkflowFeatureRequirement
inputs:
    reference:
        type:
            - string
            - File
        doc: "The reference that was used for the already-completed alignments"
    crams_to_downsample:
        type:
            type: array
            items:
                type: record
                name: crams
                fields:
                    cram:
                        type: File
                    downsample_ratio:
                        type: float
                        doc: 'the downsample ratio to use when reprocessing this CRAM'
                    contamination:
                        type: float
                        doc: 'contamination score to pass to HaplotypeCaller'
    downsample_strategy:
        type:
            - "null"
            - type: enum
              symbols: ["HighAccuracy", "ConstantMemory", "Chained"]
    downsample_seed:
        type: int?
    emit_reference_confidence:
        type:
            type: enum
            symbols: ['NONE', 'BP_RESOLUTION', 'GVCF']
    max_alternate_alleles:
        type: int?
    ploidy:
        type: int?
    read_filter:
        type: string?
    intervals:
        type:
            type: array
            items:
                type: array
                items: string
        doc: "arrays of intervals to use in each individual run of the haplotypecaller"
    qc_minimum_mapping_quality:
        type: int
    qc_minimum_base_quality:
        type: int
outputs:
    gvcfs:
        type:
            type: array
            items:
                type: array
                items: File
        outputSource: haplotype_caller/gvcf
    wgs_metrics:
        type: File[]
        outputSource: collect_wgs_metrics/wgs_metrics
steps:
    downsample:
        run: ../tools/downsample.cwl
        scatter: [sam, probability]
        scatterMethod: dotproduct
        in:
            sam:
                source: crams_to_downsample
                valueFrom: $(self.cram)
            probability:
                source: crams_to_downsample
                valueFrom: $(self.downsample_ratio)
            reference: reference
            random_seed: downsample_seed
            strategy: downsample_strategy
        out: [downsampled_sam]
    haplotype_caller:
        run: ../subworkflows/gatk_haplotypecaller_iterator.cwl
        scatter: [bam, contamination_fraction, output_prefix]
        scatterMethod: dotproduct
        in:
            reference: reference
            bam: downsample/downsampled_sam
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands:
                default: []
            intervals: intervals
            contamination_fraction:
                source: crams_to_downsample
                valueFrom: $(self.contamination)
            max_alternate_alleles: max_alternate_alleles
            ploidy: ploidy
            read_filter: read_filter
            output_prefix:
                source: downsample/downsampled_sam
                valueFrom: '$(self.nameroot + ".downsampled.")'
        out: [gvcf]
    collect_wgs_metrics:
        run: ../tools/collect_wgs_metrics.cwl
        scatter: [bam, sample_name]
        scatterMethod: dotproduct
        in:
            bam: downsample/downsampled_sam
            reference: reference
            minimum_mapping_quality: qc_minimum_mapping_quality
            minimum_base_quality: qc_minimum_base_quality
            sample_name:
                source: downsample/downsampled_sam
                valueFrom: $(self.nameroot)
        out: [wgs_metrics]
