#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "exome alignment and germline variant detection"
requirements:
    - class: SubworkflowFeatureRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/vep_custom_annotation.yml
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    bam:
        type: File
        secondaryFiles: [^.bai]
    emit_reference_confidence:
        type:
            type: enum
            symbols: ['NONE', 'BP_RESOLUTION', 'GVCF']
    gvcf_gq_bands:
        type: string[]
    intervals:
        type:
            type: array
            items:
                type: array
                items: string
    contamination_fraction:
        type: string?
    vep_cache_dir:
        type:
            - string
            - Directory
    vep_ensembl_assembly:
        type: string
        doc: "genome assembly to use in vep. Examples: GRCh38 or GRCm38"
    vep_ensembl_version:
        type: string
        doc: "ensembl version - Must be present in the cache directory. Example: 95"
    vep_ensembl_species:
        type: string
        doc: "ensembl species - Must be present in the cache directory. Examples: homo_sapiens or mus_musculus"
    vep_plugins:
        type: string[]
        default: [Downstream, Wildtype]
    synonyms_file:
        type: File?
    annotate_coding_only:
        type: boolean?
    vep_custom_annotations:
        type: ../types/vep_custom_annotation.yml#vep_custom_annotation[]
        doc: "custom type, check types directory for input format"
    limit_variant_intervals:
        type: File
    variants_to_table_fields:
        type: string[]?
        default: ['CHROM','POS','ID','REF','ALT']
    variants_to_table_genotype_fields:
        type: string[]?
    vep_to_table_fields:
        type: string[]?
    final_tsv_prefix:
        type: string?
        default: 'variants'
    filter_gnomAD_maximum_population_allele_frequency:
        type: float
        default: 0.05
outputs:
    gvcf:
        type: File[]
        outputSource: haplotype_caller/gvcf
    merged_gvcf:
        type: File
        outputSource: merge_gvcfs/merged_gvcf
        secondaryFiles: [.tbi]
steps:
    haplotype_caller:
        run: gatk_haplotypecaller_iterator_GATK4.cwl
        in:
            reference: reference
            bam: bam
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            contamination_fraction: contamination_fraction
        out:
            [gvcf]
    merge_gvcfs:
        run: ../tools/picard_merge_gvcfs.cwl
        in:
            gvcfs: haplotype_caller/gvcf
        out:
            [merged_gvcf]
