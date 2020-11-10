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
    raw_vcf:
        type: File
        outputSource: merge_vcfs/merged_vcf
        secondaryFiles: [.tbi]
    final_vcf:
        type: File
        outputSource: filter_vcf/final_vcf
        secondaryFiles: [.tbi]
    filtered_vcf:
        type: File
        outputSource: filter_vcf/filtered_vcf
        secondaryFiles: [.tbi]
    vep_summary:
        type: File
        outputSource: annotate_variants/vep_summary
    final_tsv:
        type: File
        outputSource: set_final_tsv_name/replacement
    filtered_tsv:
        type: File
        outputSource: set_filtered_tsv_name/replacement
steps:
    haplotype_caller:
        run: gatk_haplotypecaller_iterator.cwl
        in:
            reference: reference
            bam: bam
            emit_reference_confidence: emit_reference_confidence
            gvcf_gq_bands: gvcf_gq_bands
            intervals: intervals
            contamination_fraction: contamination_fraction
            ploidy: ploidy
        out:
            [gvcf]
    merge_vcfs:
        run: ../tools/picard_merge_vcfs.cwl
        in:
            gvcfs: haplotype_caller/gvcf
        out:
            [merged_vcf]
    annotate_variants:
        run: ../tools/vep.cwl
        in:
            vcf: merge_vcfs/merged_vcf
            cache_dir: vep_cache_dir
            ensembl_assembly: vep_ensembl_assembly
            ensembl_version: vep_ensembl_version
            ensembl_species: vep_ensembl_species
            synonyms_file: synonyms_file
            coding_only: annotate_coding_only
            reference: reference
            custom_annotations: vep_custom_annotations
            plugins: vep_plugins
        out:
            [annotated_vcf, vep_summary]
    bgzip_annotated_vcf:
        run: ../tools/bgzip.cwl
        in:
            file: annotate_variants/annotated_vcf
        out:
            [bgzipped_file]
    index_annotated_vcf:
        run: ../tools/index_vcf.cwl
        in:
            vcf: bgzip_annotated_vcf/bgzipped_file
        out:
            [indexed_vcf]
    filter_vcf:
        run: germline_filter_vcf.cwl
        in:
            annotated_vcf: annotate_variants/annotated_vcf
            filter_gnomAD_maximum_population_allele_frequency: filter_gnomAD_maximum_population_allele_frequency
            gnomad_field_name:
               source: vep_custom_annotations
               valueFrom: |
                 ${
                    if(self){
                         for(var i=0; i<self.length; i++){
                             if(self[i].annotation.gnomad_filter){
                                 return(self[i].annotation.name + '_AF');
                             }
                         }
                     }
                     return('gnomAD_AF');
                 }
            limit_variant_intervals: limit_variant_intervals
            reference: reference
        out:
            [filtered_vcf, final_vcf]
    filtered_variants_to_table:
        run: ../tools/variants_to_table.cwl
        in:
            reference: reference
            vcf: filter_vcf/filtered_vcf
            fields: variants_to_table_fields
            genotype_fields: variants_to_table_genotype_fields
        out:
            [variants_tsv]
    filtered_add_vep_fields_to_table:
        run: ../tools/add_vep_fields_to_table.cwl
        in:
            vcf: filter_vcf/filtered_vcf
            vep_fields: vep_to_table_fields
            tsv: filtered_variants_to_table/variants_tsv
            prefix: final_tsv_prefix
        out:
            [annotated_variants_tsv]
    set_filtered_tsv_name:
        run: ../tools/staged_rename.cwl
        in:
            original: filtered_add_vep_fields_to_table/annotated_variants_tsv
            name:
                valueFrom: 'annotated.filtered.tsv'
        out:
             [replacement]
    final_variants_to_table:
        run: ../tools/variants_to_table.cwl
        in:
            reference: reference
            vcf: filter_vcf/final_vcf
            fields: variants_to_table_fields
            genotype_fields: variants_to_table_genotype_fields
        out:
            [variants_tsv]
    final_add_vep_fields_to_table:
        run: ../tools/add_vep_fields_to_table.cwl
        in:
            vcf: filter_vcf/final_vcf
            vep_fields: vep_to_table_fields
            tsv: final_variants_to_table/variants_tsv
            prefix: final_tsv_prefix
        out:
            [annotated_variants_tsv]
    set_final_tsv_name:
        run: ../tools/staged_rename.cwl
        in:
            original: final_add_vep_fields_to_table/annotated_variants_tsv
            name:
                valueFrom: 'annotated.filtered.final.tsv'
        out:
             [replacement]
