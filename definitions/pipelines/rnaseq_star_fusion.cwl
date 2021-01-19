#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "STAR-RNA-Seq alignment and transcript/gene abundance workflow"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
inputs:
    instrument_data_bams:
        type: File[]
    outsam_attrrg_line:
        type: string[]
    star_genome_dir:
        type: Directory
    star_fusion_genome_dir:
        type: Directory
    cdna_fasta:
        type: File
    reference_fasta:
        type: File
        secondaryFiles: [.fai, ^.dict]
    gtf_file:
        type: File
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
        type:
          - "null"
          - type: enum
            symbols: ["first", "second", "unstranded"]
    refFlat:
        type: File
    ribosomal_intervals:
        type: File
    sample_name:
        type: string
outputs:
    cram:
        type: File
        outputSource: index_cram/indexed_cram
        secondaryFiles: [.crai, ^.crai]
    star_fusion_out:
        type: File
        outputSource: star_align_fusion/chim_junc
    star_junction_out:
        type: File
        outputSource: star_align_fusion/splice_junction_out
    star_fusion_log:
        type: File
        outputSource: star_align_fusion/log_final
    star_fusion_predict:
        type: File
        outputSource: star_fusion_detect/fusion_predictions
    star_fusion_abridge:
        type: File
        outputSource: star_fusion_detect/fusion_abridged
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
    strand_info:
        type: File[]
        outputSource: strandedness_check/strandedness_check
    bamcoverage_bigwig:
        type: File
        outputSource: cgpbigwig_bamcoverage/outfile
steps:
    bam_to_trimmed_fastq:
        run: ../subworkflows/bam_to_trimmed_fastq.cwl
        scatter: [bam]
        scatterMethod: dotproduct
        in:
            bam: instrument_data_bams
            adapters: trimming_adapters
            adapter_trim_end: trimming_adapter_trim_end
            adapter_min_overlap: trimming_adapter_min_overlap
            max_uncalled: trimming_max_uncalled
            min_readlength: trimming_min_readlength
        out:
            [fastqs, fastq1, fastq2]
    strandedness_check:
        run: ../tools/strandedness_check.cwl
        scatter: [reads1, reads2]
        scatterMethod: dotproduct
        in:
            gtf_file: gtf_file
            kallisto_index: kallisto_index
            cdna_fasta: cdna_fasta
            reads1: bam_to_trimmed_fastq/fastq1
            reads2: bam_to_trimmed_fastq/fastq2
        out:
            [strandedness_check]
    star_align_fusion:
        run: ../tools/star_align_fusion.cwl
        in:
            outsam_attrrg_line: outsam_attrrg_line
            star_genome_dir: star_genome_dir
            gtf_file: gtf_file
            fastq:
                source: bam_to_trimmed_fastq/fastq1
                linkMerge: merge_flattened
            fastq2:
                source: bam_to_trimmed_fastq/fastq2
                linkMerge: merge_flattened
        out:
            [aligned_bam, chim_junc, splice_junction_out,log_final]
    star_fusion_detect:
        run: ../tools/star_fusion_detect.cwl
        in:
            star_fusion_genome_dir: star_fusion_genome_dir
            junction_file: star_align_fusion/chim_junc
        out:
            [fusion_predictions,fusion_abridged]
    kallisto:
        run: ../tools/kallisto.cwl
        in:
            kallisto_index: kallisto_index
            strand: strand
            fastqs: bam_to_trimmed_fastq/fastqs
        out:
            [expression_transcript_table,expression_transcript_h5,fusion_evidence]
    transcript_to_gene:
        run: ../tools/transcript_to_gene.cwl
        in:
            transcript_table_h5: kallisto/expression_transcript_h5
            gene_transcript_lookup_table: gene_transcript_lookup_table
        out:
            [gene_abundance]
    sort_bam:
        run: ../tools/samtools_sort.cwl
        in:
            input_bam: star_align_fusion/aligned_bam
        out:
            [sorted_bam]
    mark_dup:
        run: ../tools/mark_duplicates_and_sort.cwl
        in:
            bam: sort_bam/sorted_bam
            input_sort_order:
                default: "coordinate"
        out:
            [sorted_bam, metrics_file]
    index_bam:
        run: ../tools/index_bam.cwl
        in:
            bam: mark_dup/sorted_bam
        out:
            [indexed_bam]
    stringtie:
        run: ../tools/stringtie.cwl
        in:
            bam: mark_dup/sorted_bam
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
            bam: mark_dup/sorted_bam
        out:
            [metrics, chart]
    bam_to_cram:
        run: ../tools/bam_to_cram.cwl
        in:
          reference: reference_fasta
          bam: index_bam/indexed_bam
        out:
            [cram]
    index_cram:
        run: ../tools/index_cram.cwl
        in:
            cram: bam_to_cram/cram
        out:
            [indexed_cram]
    cgpbigwig_bamcoverage:
        run: ../tools/bam_to_bigwig.cwl
        in:
            bam: index_bam/indexed_bam
            reference: reference_fasta
        out:
            [outfile]
