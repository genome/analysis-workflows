#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "RNA-Seq alignment and transcript/gene abundance workflow"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
    - class: InlineJavascriptRequirement
inputs:
    instrument_data_bams:
        type: File[]
    outSAMattrRGline:
        type:string[]
    stargenomeDir:
        type: Directory
    gtf_file:
        type: File
    sample_name:
        type: string
    outFileNamePrefix:
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
    strand:
       type: string?
    refFlat:
        type: File
    ribosomal_intervals:
        type: File
outputs:
    final_bam:
        type: File
        outputSource: index_bam/indexed_bam
        secondaryFiles: [.bai]
    stringtie_transcript_gtf:
        type: File
        outputSource: stringtie/transcript_gtf
    stringtie_gene_expression_tsv:
        type: File
        outputSource: stringtie/gene_expression_tsv
    transcript_abundance_tsv:
        type: File
        outputSource: kallisto/expression_transcript_table
    transcript_abundance_h5:
        type: File
        outputSource: kallisto/expression_transcript_h5
    gene_abundance:
        type: File
        outputSource: transcript_to_gene/gene_abundance
    metrics:
        type: File
        outputSource: generate_qc_metrics/metrics
    chart:
        type: File
        outputSource: generate_qc_metrics/chart
    fusion_evidence:
        type: File
        outputSource: kallisto/fusion_evidence
steps:
    bam_to_trimmed_fastq_and_star_alignments:
        run: ../subworkflows/bam_to_trimmed_fastq_and_star_alignments.cwl
        scatter: [bam, outSAMattrRGline]
        scatterMethod: dotproduct
        in:
            bam: instrument_data_bams
            outSAMattrRGline: outSAMattrRGline
            outFileNamePrefix: sample_name
            stargenomeDir: stargenomeDir
            adapters: trimming_adapters
            adapter_trim_end: trimming_adapter_trim_end
            adapter_min_overlap: trimming_adapter_min_overlap
            max_uncalled: trimming_max_uncalled
            min_readlength: trimming_min_readlength
        out:
            [fastqs,aligned_bam]
    kallisto:
        run: ../tools/kallisto.cwl
        in:
            kallisto_index: kallisto_index
            strand: strand
            fastqs:
                source: bam_to_trimmed_fastq_and_star_alignments/fastqs
                valueFrom: |
                    ${
                      for(var i=0;i<self.length;i++){self[i] = self[i].reverse()}
                      return(self)
                     }
        out:
            [expression_transcript_table,expression_transcript_h5,fusion_evidence]
    transcript_to_gene:
        run: ../tools/transcript_to_gene.cwl
        in:
            transcript_table_h5: kallisto/expression_transcript_h5
            gene_transcript_lookup_table: gene_transcript_lookup_table
        out:
            [gene_abundance]
    merge:
        run: ../tools/merge_bams.cwl
        in:
            bams: bam_to_trimmed_fastq_and_star_alignments/aligned_bam
        out:
            [merged_bam]
    index_bam:
        run: ../tools/index_bam.cwl
        in:
            bam: merge/merged_bam
        out:
            [indexed_bam]
    stringtie:
        run: ../tools/stringtie.cwl
        in:
            bam: merge/merged_bam
            reference_annotation: gtf_file
            sample_name: sample_name
            strand: strand
        out:
            [transcript_gtf,gene_expression_tsv]
    generate_qc_metrics:
        run: ../tools/generate_qc_metrics.cwl
        in:
            refFlat: refFlat
            ribosomal_intervals: ribosomal_intervals
            strand: strand
            bam: index_bam/indexed_bam
        out:
            [metrics, chart]
