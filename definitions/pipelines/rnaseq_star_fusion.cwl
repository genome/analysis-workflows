#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "STAR-RNA-Seq alignment and transcript/gene abundance workflow"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
inputs:
    unaligned:
        type: ../types/sequence_data.yml#sequence_data[]
        doc: |
            Raw data from rna sequencing; this custom type holds both the data
            file(s) and readgroup information. Data file(s) may be either a bam
            file, or paired fastqs. Readgroup information should be given as a
            series of key:value pairs, each separated by a space. This means
            that spaces within a value must be double quoted. The first key
            must be ID; consult the read group description in the header
            section of the SAM file specification for other, optional keys.
            Below is an example of an element of the input array:
            readgroup: "ID:xxx PU:xxx SM:xxx LB:xxx PL:ILLUMINA CN:WUGSC"
            sequence:
                fastq1:
                    class: File
                    path: /path/to/reads1.fastq
                fastq2:
                    class: File
                    path: /path/to/reads2.fastq
                OR
                bam:
                    class: File
                    path: /path/to/reads.bam
    star_genome_dir:
        type: Directory
    star_fusion_genome_dir:
        type: Directory
    cdna_fasta:
        type: File
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    reference_annotation:
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
    unzip_fastqs:
        type: boolean?
        default: true
    examine_coding_effect:
        type: boolean?
    fusioninspector_mode:
        type:
          - "null"
          - type: enum
            symbols: ["inspect", "validate"]
    agfusion_database:
        type: File
    agfusion_annotate_noncanonical:
        type: boolean?

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
    final_bam:
        type: File
        outputSource: index_bam/indexed_bam
        secondaryFiles: [.bai]
    annotated_fusion_predictions:
        type: Directory
        outputSource: agfusion/annotated_fusion_predictions
    coding_region_effects:
        type: File?
        outputSource: star_fusion_detect/coding_region_effects
    fusioninspector_evidence:
        type: File[]?
        outputSource: star_fusion_detect/fusioninspector_evidence
steps:
    sequence_to_trimmed_fastq:
        
        scatter: [unaligned]
        scatterMethod: dotproduct
        run: ../subworkflows/sequence_to_trimmed_fastq.cwl
        in:
            unaligned: unaligned
            adapters: trimming_adapters
            adapter_trim_end: trimming_adapter_trim_end
            adapter_min_overlap: trimming_adapter_min_overlap
            max_uncalled: trimming_max_uncalled
            min_readlength: trimming_min_readlength
            unzip_fastqs: unzip_fastqs
        out:
            [fastqs, fastq1, fastq2]
    strandedness_check:
        run: ../tools/strandedness_check.cwl
        scatter: [reads1, reads2]
        scatterMethod: dotproduct
        in:
            reference_annotation: reference_annotation
            kallisto_index: kallisto_index
            cdna_fasta: cdna_fasta
            reads1: sequence_to_trimmed_fastq/fastq1
            reads2: sequence_to_trimmed_fastq/fastq2
        out:
            [strandedness_check]
    star_align_fusion:
        run: ../tools/star_align_fusion.cwl
        in:
            outsam_attrrg_line:
                source: unaligned
                valueFrom: $(self.map(seqtype => seqtype.readgroup))
            star_genome_dir: star_genome_dir
            reference_annotation: reference_annotation
            fastq:
                source: sequence_to_trimmed_fastq/fastq1
                linkMerge: merge_flattened
            fastq2:
                source: sequence_to_trimmed_fastq/fastq2
                linkMerge: merge_flattened
        out:
            [aligned_bam, chim_junc, splice_junction_out,log_final]
    star_fusion_detect:
        run: ../tools/star_fusion_detect.cwl
        in:
            star_fusion_genome_dir: star_fusion_genome_dir
            junction_file: star_align_fusion/chim_junc
            examine_coding_effect: examine_coding_effect
            fusioninspector_mode: fusioninspector_mode
            fastq:
                source: sequence_to_trimmed_fastq/fastq1
                linkMerge: merge_flattened
            fastq2:
                source: sequence_to_trimmed_fastq/fastq2
                linkMerge: merge_flattened
        out:
            [fusion_predictions,fusion_abridged, coding_region_effects, fusioninspector_evidence]
    kallisto:
        run: ../tools/kallisto.cwl
        in:
            kallisto_index: kallisto_index
            strand: strand
            fastqs: sequence_to_trimmed_fastq/fastqs
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
            reference_annotation: reference_annotation
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
          reference: reference
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
            bam: mark_dup/sorted_bam
            reference: reference
        out:
            [outfile]
    agfusion:
        run: ../tools/agfusion.cwl
        in:
            fusion_predictions: star_fusion_detect/fusion_predictions
            agfusion_database: agfusion_database
            annotate_noncanonical: agfusion_annotate_noncanonical
        out:
            [annotated_fusion_predictions]
