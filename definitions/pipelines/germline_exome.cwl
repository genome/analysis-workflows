#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "exome alignment and germline variant detection"
requirements:
    - class: SubworkflowFeatureRequirement
inputs:
    reference: string
    bams:
        type: File[]
    readgroups:
        type: string[]
    mills:
        type: File
        secondaryFiles: [.tbi]
    known_indels:
        type: File
        secondaryFiles: [.tbi]
    dbsnp_vcf:
        type: File
        secondaryFiles: [.tbi]
    bqsr_intervals:
        type: string[]?
    bait_intervals:
        type: File
    target_intervals:
        type: File
    per_target_intervals:
        type: File
    per_target_bait_intervals:
        type: File
    per_base_intervals:
        type: File
    per_base_bait_intervals:
        type: File
    omni_vcf:
        type: File
        secondaryFiles: [.tbi]
    picard_metric_accumulation_level:
        type: string
    emit_reference_confidence:
        type: string
    gvcf_gq_bands:
        type: string[]
    intervals:
        type:
            type: array
            items:
                type: array
                items: string
    vep_cache_dir:
        type: string?
    synonyms_file:
        type: File?
    annotate_coding_only:
        type: boolean?
    hgvs_annotation:
        type: boolean?
    custom_gnomad_vcf:
        type: File?
        secondaryFiles: [.tbi]
    qc_minimum_mapping_quality:
        type: int?
    qc_minimum_base_quality:
        type: int?
    output_dir:
        type: string
outputs:
    cram:
        type: File
        outputSource: alignment_and_qc/cram
    mark_duplicates_metrics:
        type: File
        outputSource: alignment_and_qc/mark_duplicates_metrics
    insert_size_metrics:
        type: File
        outputSource: alignment_and_qc/insert_size_metrics
    insert_size_histogram:
        type: File
        outputSource: alignment_and_qc/insert_size_histogram
    alignment_summary_metrics:
        type: File
        outputSource: alignment_and_qc/alignment_summary_metrics
    hs_metrics:
        type: File
        outputSource: alignment_and_qc/hs_metrics
    per_target_coverage_metrics:
        type: File?
        outputSource: alignment_and_qc/per_target_coverage_metrics
    per_target_hs_metrics:
        type: File?
        outputSource: alignment_and_qc/per_target_hs_metrics
    per_base_coverage_metrics:
        type: File?
        outputSource: alignment_and_qc/per_base_coverage_metrics
    per_base_hs_metrics:
        type: File?
        outputSource: alignment_and_qc/per_base_hs_metrics
    flagstats:
        type: File
        outputSource: alignment_and_qc/flagstats
    verify_bam_id_metrics:
        type: File
        outputSource: alignment_and_qc/verify_bam_id_metrics
    verify_bam_id_depth:
        type: File
        outputSource: alignment_and_qc/verify_bam_id_depth
    gvcf:
        type: File[]
        outputSource: detect_variants/gvcf
    final_vcf:
        type: File
        outputSource: detect_variants/final_vcf
        secondaryFiles: [.tbi]
    coding_vcf:
        type: File
        outputSource: detect_variants/coding_vcf
        secondaryFiles: [.tbi]
    limited_vcf:
        type: File
        outputSource: detect_variants/limited_vcf
        secondaryFiles: [.tbi]
    vep_summary:
        type: File
        outputSource: detect_variants/vep_summary
    final_outputs:
        type: string[]
        outputSource: gatherer/gathered_files
steps:
    alignment_and_qc:
        run: exome_alignment.cwl
        in:
            reference: reference
            bams: bams
            readgroups: readgroups
            mills: mills
            known_indels: known_indels
            dbsnp_vcf: dbsnp_vcf
            bqsr_intervals: bqsr_intervals
            bait_intervals: bait_intervals
            target_intervals: target_intervals
            per_target_intervals: per_target_intervals
            per_target_bait_intervals: per_target_bait_intervals
            per_base_intervals: per_base_intervals
            per_base_bait_intervals: per_base_bait_intervals
            omni_vcf: omni_vcf
            picard_metric_accumulation_level: picard_metric_accumulation_level   
            minimum_mapping_quality: qc_minimum_mapping_quality
            minimum_base_quality: qc_minimum_base_quality
        out:
            [cram, mark_duplicates_metrics, insert_size_metrics, insert_size_histogram, alignment_summary_metrics, hs_metrics, per_target_coverage_metrics, per_target_hs_metrics, per_base_coverage_metrics, per_base_hs_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
    extract_freemix:
        in:
            verify_bam_id_metrics: alignment_and_qc/verify_bam_id_metrics
        out:
            [freemix_score]
        run:
            class: ExpressionTool
            requirements:
                - class: InlineJavascriptRequirement
            inputs:
                verify_bam_id_metrics:
                    type: File
                    inputBinding:
                        loadContents: true
            outputs:
                freemix_score:
                    type: string?
            expression: |
                        ${
                            var metrics = inputs.verify_bam_id_metrics.contents.split("\n");
                            if ( metrics[0].split("\t")[6] == 'FREEMIX' ) {
                                return {'freemix_score': metrics[1].split("\t")[6]};
                            } else {
                                return {'freemix_score:': null };
                            }
                        }
    detect_variants:
        run: ../subworkflows/germline_detect_variants.cwl
        in:
            reference: reference
            cram: alignment_and_qc/cram
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            contamination_fraction: extract_freemix/freemix_score
            vep_cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            hgvs: hgvs_annotation
            custom_gnomad_vcf: custom_gnomad_vcf
            limit_variant_intervals: per_base_intervals
        out:
            [gvcf, final_vcf, coding_vcf, limited_vcf, vep_summary]
    gatherer:
        run: ../tools/gatherer.cwl
        in:
            output_dir: output_dir
            all_files:
                source: [alignment_and_qc/per_target_hs_metrics, alignment_and_qc/hs_metrics, alignment_and_qc/per_base_hs_metrics, alignment_and_qc/alignment_summary_metrics, detect_variants/coding_vcf, detect_variants/gvcf, detect_variants/limited_vcf, detect_variants/vep_summary, detect_variants/final_vcf, alignment_and_qc/flagstats, alignment_and_qc/insert_size_metrics, alignment_and_qc/mark_duplicates_metrics, alignment_and_qc/per_base_coverage_metrics, alignment_and_qc/verify_bam_id_metrics, alignment_and_qc/verify_bam_id_depth, alignment_and_qc/insert_size_histogram, alignment_and_qc/per_target_coverage_metrics, alignment_and_qc/cram]
                valueFrom: ${
                                function flatten(inArr, outArr) {
                                    var arrLen = inArr.length;
                                    for (var i = 0; i < arrLen; i++) {
                                        if (Array.isArray(inArr[i])) {
                                            flatten(inArr[i], outArr);
                                        }
                                        else {
                                            outArr.push(inArr[i]);
                                        }
                                    }
                                    return outArr;
                                }
                                var no_secondaries = flatten(self, []);
                                var all_files = []; 
                                var arrLen = no_secondaries.length;
                                for (var i = 0; i < arrLen; i++) {
                                    all_files.push(no_secondaries[i]);
                                    var secondaryLen = no_secondaries[i].secondaryFiles.length;
                                    for (var j = 0; j < secondaryLen; j++) {
                                        all_files.push(no_secondaries[i].secondaryFiles[j]);
                                    }
                                }
                                return all_files;
                            }
        out: [gathered_files]
