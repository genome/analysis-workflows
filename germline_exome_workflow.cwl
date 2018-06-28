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
    dbsnp:
        type: File
        secondaryFiles: [.tbi]
    bqsr_intervals:
        type: string[]?
    bait_intervals:
        type: File
    target_intervals:
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
    coding_only:
        type: boolean?
    hgvs_annotation:
        type: boolean?
    custom_gnomad_vcf:
        type: File?
        secondaryFiles: [.tbi]
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
    alignment_summary_metrics:
        type: File
        outputSource: alignment_and_qc/alignment_summary_metrics
    hs_metrics:
        type: File
        outputSource: alignment_and_qc/hs_metrics
    per_target_coverage_metrics:
        type: File?
        outputSource: alignment_and_qc/per_target_coverage_metrics
    per_base_coverage_metrics:
        type: File?
        outputSource: alignment_and_qc/per_base_coverage_metrics
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
        outputSource: haplotype_caller/gvcf
    final_vcf:
        type: File
        outputSource: index_annotated_vcf/indexed_vcf
        secondaryFiles: [.tbi]
    vep_summary:
        type: File
        outputSource: annotate_variants/vep_summary
steps:
    alignment_and_qc:
        run: exome_alignment.cwl
        in:
            reference: reference
            bams: bams
            readgroups: readgroups
            mills: mills
            known_indels: known_indels
            dbsnp: dbsnp
            bqsr_intervals: bqsr_intervals
            bait_intervals: bait_intervals
            target_intervals: target_intervals
            omni_vcf: omni_vcf
            picard_metric_accumulation_level: picard_metric_accumulation_level   
        out:
            [cram, mark_duplicates_metrics, insert_size_metrics, alignment_summary_metrics, hs_metrics, per_target_coverage_metrics, per_base_coverage_metrics, flagstats, verify_bam_id_metrics, verify_bam_id_depth]
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
    haplotype_caller:
        run: detect_variants/gatk_haplotypecaller_iterator.cwl
        in:
            reference: reference
            cram: alignment_and_qc/cram
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            contamination_fraction: extract_freemix/freemix_score
        out:
            [gvcf]
    genotype_gvcfs:
        run: detect_variants/gatk_genotypegvcfs.cwl
        in:
            reference: reference
            gvcfs: haplotype_caller/gvcf
        out:
            [genotype_vcf]
    annotate_variants:
        run: detect_variants/vep.cwl
        in:
            vcf: genotype_gvcfs/genotype_vcf
            cache_dir: vep_cache_dir
            synonyms_file: synonyms_file
            coding_only: coding_only
            hgvs: hgvs_annotation
            reference: reference
            custom_gnomad_vcf: custom_gnomad_vcf
        out:
            [annotated_vcf, vep_summary]
    bgzip_annotated_vcf:
        run: detect_variants/bgzip.cwl
        in:
            file: annotate_variants/annotated_vcf
        out:
            [bgzipped_file]
    index_annotated_vcf:
        run: detect_variants/index.cwl
        in:
            vcf: bgzip_annotated_vcf/bgzipped_file
        out:
            [indexed_vcf]
