#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "RNA-Seq alignment and transcript/gene abundance workflow - first-stranded data"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
    - class: InlineJavascriptRequirement
inputs:
    reference_index:
        type: File #this requires an extra file with the basename
        secondaryFiles: [".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2", ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"]
    reference_annotation:
        type: File
    instrument_data_bams:
        type: File[]
    read_group_id:
        type: string[]
    read_group_fields:
        type:
            type: array
            items:
                type: array
                items: string
    sample_name:
        type: string
    trimming_adapters:
        type: File
    trimming_adapter_trim_end:
        type: string
    trimming_adapter_min_overlap:
        type: int
    trimming_max_uncalled:
        type: int
    trimming_min_readlength:
        type: int
    kallisto_index: 
       type: File
    gene_transcript_lookup_table:
       type: File
    firststrand:
       type: boolean?
    secondstrand:
       type: boolean?
outputs:
    final_bam:
        type: File
        outputSource: index_bam/indexed_bam
        secondaryFiles: [.bai, ^.bai]
    gtf:
        type: File
        outputSource: stringtie/gtf
    transcript_abundance_tsv:
        type: File
        outputSource: kallisto/expression_transcript_table
    transcript_abundance_h5:
        type: File
        outputSource: kallisto/expression_transcript_h5
    gene_abundance:
        type: File
        outputSource: transcript_to_gene/gene_abundance
steps:
    bam_to_trimmed_fastq_and_hisat_alignments:
        run: bam_to_trimmed_fastq_and_hisat_alignments.cwl
        scatter: [bam, read_group_id, read_group_fields]
        scatterMethod: dotproduct
        in:
            bam: instrument_data_bams
            read_group_id: read_group_id
            read_group_fields: read_group_fields
            adapters: trimming_adapters
            adapter_trim_end: trimming_adapter_trim_end
            adapter_min_overlap: trimming_adapter_min_overlap
            max_uncalled: trimming_max_uncalled
            min_readlength: trimming_min_readlength
            reference_index: reference_index
            firststrand: firststrand
            secondstrand: secondstrand
        out:
            [fastqs,aligned_bam]
    kallisto:
        run: kallisto.cwl
        in:
            kallisto_index: kallisto_index
            firststrand: firststrand
            secondstrand: secondstrand
            fastqs: 
                source: bam_to_trimmed_fastq_and_hisat_alignments/fastqs
                valueFrom: |
                    ${
                      for(var i=0;i<self.length;i++){self[i] = self[i].reverse()}
                      return(self)                      
                     }
        out:
            [expression_transcript_table,expression_transcript_h5]
    transcript_to_gene:
        run: transcript_to_gene.cwl
        in:
            transcript_table_h5: kallisto/expression_transcript_h5
            gene_transcript_lookup_table: gene_transcript_lookup_table
        out:
            [gene_abundance]
    merge:
        run: merge.cwl
        in:
            bams: bam_to_trimmed_fastq_and_hisat_alignments/aligned_bam
        out:
            [merged_bam]
    index_bam:
        run: ../detect_variants/index_bam.cwl
        in:
            bam: merge/merged_bam
        out:
            [indexed_bam]
    stringtie:
        run: stringtie.cwl
        in:
            bam: merge/merged_bam
            reference_annotation: reference_annotation
            sample_name: sample_name
        out:
            [gtf]
