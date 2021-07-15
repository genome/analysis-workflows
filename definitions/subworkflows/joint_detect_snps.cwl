#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "joint germline snp variant detection"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/vep_custom_annotation.yml
    - class: StepInputExpressionRequirement
    - class: InlineJavascriptRequirement
    - class: ScatterFeatureRequirement
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    bams:
        type: File[]
        secondaryFiles: [^.bai]
    sample_names:
        type: string[]
    gvcf_gq_bands:
        type: string[]
    intervals:
        type:
            type: array
            items:
                type: array
                items: string
    contamination_fraction:
        type: string[]
    ploidy:
        type: int?
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
        default: [Frameshift, Wildtype]
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
        type: string[]
        default: ['CHROM','POS','ID','REF','ALT']
    variants_to_table_genotype_fields:
        type: string[]
    vep_to_table_fields:
        type: string[]
    final_tsv_prefix:
        type: string?
        default: 'variants'
    gnomad_max_pop_af:
        type: float
        default: 0.05
    min_conf_call:
        type: float?
outputs:
    sample_gvcfs:
        type: File[]
        outputSource: per_sample_merge_gvcfs/gvcf
    raw_vcf:
        type: File
        outputSource: genotype/raw_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: genotype/final_vcf
        secondaryFiles: [.tbi]
    filtered_vcf:
        type: File
        outputSource: genotype/filtered_vcf
        secondaryFiles: [.tbi]
    vep_summary:
        type: File
        outputSource: genotype/vep_summary
    final_tsv:
        type: File
        outputSource: genotype/final_tsv
    filtered_tsv:
        type: File
        outputSource: genotype/filtered_tsv
    all_staged:
        type: Directory
        outputSource: stage_all/gathered_directory
steps:
    per_sample_make_gvcfs:
        scatter: [bam, contamination_fraction]
        scatterMethod: dotproduct
        run: gatk_haplotypecaller_iterator.cwl
        in:
            reference: reference
            bam: bams
            emit_reference_confidence:
                 default: 'GVCF'
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            contamination_fraction: contamination_fraction
            ploidy: ploidy
        out:
            [gvcf]
    per_sample_merge_gvcfs:
        scatter: [gvcfs, output_file_name]
        scatterMethod: dotproduct
        run: ../tools/combine_gvcfs.cwl
        in:
            reference: reference
            gvcfs: per_sample_make_gvcfs/gvcf
            output_file_name:
                source: [sample_names]
                valueFrom: "$(self).merged.g.vcf.gz"
        out:
            [gvcf]
    genotype:
        run: joint_genotype.cwl
        in:
            reference: reference
            gvcfs:
                source: [per_sample_merge_gvcfs/gvcf]
                linkMerge: merge_flattened
            intervals: intervals
            vep_cache_dir: vep_cache_dir
            vep_ensembl_assembly: vep_ensembl_assembly
            vep_ensembl_version: vep_ensembl_version
            vep_ensembl_species: vep_ensembl_species
            vep_plugins: vep_plugins
            synonyms_file: synonyms_file
            annotate_coding_only: annotate_coding_only
            vep_custom_annotations: vep_custom_annotations
            roi_intervals: limit_variant_intervals
            variants_to_table_fields: variants_to_table_fields
            variants_to_table_genotype_fields: variants_to_table_genotype_fields
            vep_to_table_fields: vep_to_table_fields
            final_tsv_prefix: final_tsv_prefix
            gnomad_max_pop_af: gnomad_max_pop_af
            min_conf_call: min_conf_call
        out:
            [raw_vcf, annotated_vcf, final_vcf, filtered_vcf, vep_summary, final_tsv, filtered_tsv]
    stage_gvcf:
        run: ../tools/gather_to_sub_directory.cwl
        in:
            outdir:
                valueFrom: "gvcfs"
            files:
                source: [per_sample_merge_gvcfs/gvcf]
                linkMerge: merge_flattened
        out:
            [gathered_directory]

    stage_gatk:
        run: ../tools/gather_to_sub_directory.cwl
        in:
            outdir:
                valueFrom: "gatk"
            files:
                source: [genotype/raw_vcf, genotype/annotated_vcf, genotype/final_vcf, genotype/filtered_vcf, genotype/vep_summary, genotype/final_tsv, genotype/filtered_tsv]
                linkMerge: merge_flattened
            directory: stage_gvcf/gathered_directory
        out:
            [gathered_directory]
    stage_all:
        run: ../tools/gather_to_sub_directory_dirs.cwl
        in:
             outdir:
                 valueFrom: "SNP_pipeline"
             directories:
                 source: [stage_gatk/gathered_directory]
                 linkMerge: merge_flattened
        out:
            [gathered_directory]
