#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Gathered Downsample and HaplotypeCaller"
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
    variant_index_type:
        type:
            - 'null'
            - type: enum
              symbols: ['DYNAMIC_SEEK', 'DYNAMIC_SIZE', 'LINEAR', 'INTERVAL']
    variant_index_parameter:
        type: string?
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
    result_directories:
        type: Directory[]
        outputSource: gather_results/gathered_directory
steps:
    downsample_and_recall:
        run: downsample_and_recall.cwl
        in:
            reference: reference
            crams_to_downsample: crams_to_downsample
            downsample_strategy: downsample_strategy
            downsample_seed: downsample_seed
            emit_reference_confidence: emit_reference_confidence
            max_alternate_alleles: max_alternate_alleles
            variant_index_type: variant_index_type
            variant_index_parameter: variant_index_parameter
            read_filter: read_filter
            intervals: intervals
            qc_minimum_mapping_quality: qc_minimum_mapping_quality
            qc_minimum_base_quality: qc_minimum_base_quality
        out: [gvcfs, wgs_metrics]
    join_results:
        in:
            gvcfs: downsample_and_recall/gvcfs
            wgs_metrics: downsample_and_recall/wgs_metrics
        out:
            [results]
        scatter: [gvcfs, wgs_metrics]
        scatterMethod: dotproduct
        run:
            class: ExpressionTool
            requirements:
                - class: InlineJavascriptRequirement
            inputs:
                gvcfs:
                    type: File[]
                wgs_metrics:
                    type: File
            outputs:
                results:
                    type: File[]
            expression: |
                ${
                    var results = [inputs.wgs_metrics];
                    results = results.concat(inputs.gvcfs);
                    return {'results': results};
                }
    gather_results:
        run: ../tools/gather_to_sub_directory.cwl
        scatter: [outdir, files]
        scatterMethod: dotproduct
        in:
            outdir:
                source: crams_to_downsample
                valueFrom: $(self.cram.nameroot)
            files: join_results/results
        out: [gathered_directory]
