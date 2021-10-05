#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "STAR-RNA-Seq alignment and transcript/gene abundance workflow with Xenosplit"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: ScatterFeatureRequirement
    - class: StepInputExpressionRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml
inputs:
    reference:
        type:
            - string
            - File
        secondaryFiles: [.fai, ^.dict]
    instrument_data_bams:
        type: ../types/sequence_data.yml#sequence_data[]
    outsam_attrrg_line:
        type: string[]
    graft_star_genome_dir:
        type: Directory
        doc: "Location of the STAR reference for the graft species (usually human)"
    graft_outfile_name_prefix:
        type: string?
        default: "Graft_"
        doc: "Prefix name for the graft bam (often 'human')"
    host_star_genome_dir:
        type: Directory
        doc: "Location of the STAR reference for the host species (usually mouse)"
    host_outfile_name_prefix:
        type: string?
        default: "Host_"
        doc: "Prefix name for the host bam (often 'mouse')"
    star_fusion_genome_dir:
        type: Directory
        doc: "star fusion directory for the graft species (often human) - fusions are only called on the xenograft-filtered bam"
    graft_gtf_file:
        type: File
        doc: "GTF file corresponding the to the graft species (often human)"
    host_gtf_file:
        type: File
        doc: "GTF file corresponding the to the host species (often mouse)"
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
    sample_name:
        type: string
outputs:
    final_bam:
        type: File
        outputSource: index_bam/indexed_bam
        secondaryFiles: [.bai]
    star_fusion_out:
        type: File
        outputSource: graft_star_align_fusion/chim_junc
    star_junction_out:
        type: File
        outputSource: graft_star_align_fusion/splice_junction_out
    star_fusion_log:
        type: File
        outputSource: graft_star_align_fusion/log_final
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
    xenosplit_statistics:
        type: File
        outputSource: xenosplit/xenosplit_statistics
    bamcoverage_bigwig:
        type: File
        outputSource: cgpbigwig_bamcoverage/outfile
steps:
    bam_to_trimmed_fastq:
        run: ../subworkflows/bam_to_trimmed_fastq.cwl
        scatter: [unaligned]
        scatterMethod: dotproduct
        in:
            unaligned: instrument_data_bams
            adapters: trimming_adapters
            adapter_trim_end: trimming_adapter_trim_end
            adapter_min_overlap: trimming_adapter_min_overlap
            max_uncalled: trimming_max_uncalled
            min_readlength: trimming_min_readlength
        out:
            [fastqs, fastq1, fastq2]
    graft_star_align_fusion:
        run: ../tools/star_align_fusion.cwl
        in:
            outsam_attrrg_line: outsam_attrrg_line
            star_genome_dir: graft_star_genome_dir
            outfile_name_prefix: graft_outfile_name_prefix
            gtf_file: graft_gtf_file
            fastq:
                source: bam_to_trimmed_fastq/fastq1
                linkMerge: merge_flattened
            fastq2:
                source: bam_to_trimmed_fastq/fastq2
                linkMerge: merge_flattened
        out:
            [aligned_bam, chim_junc, splice_junction_out,log_final]
    host_star_align_fusion:
        run: ../tools/star_align_fusion.cwl
        in:
            outsam_attrrg_line: outsam_attrrg_line
            star_genome_dir: host_star_genome_dir
            outfile_name_prefix: host_outfile_name_prefix
            gtf_file: host_gtf_file
            fastq:
                source: bam_to_trimmed_fastq/fastq1
                linkMerge: merge_flattened
            fastq2:
                source: bam_to_trimmed_fastq/fastq2
                linkMerge: merge_flattened
        out:
            [aligned_bam, chim_junc, splice_junction_out,log_final]   
    xenosplit:
        run: ../tools/xenosplit.cwl
        in:
            graftbam: graft_star_align_fusion/aligned_bam
            hostbam: host_star_align_fusion/aligned_bam
        out:
            [graft_bam, xenosplit_statistics]
    graftbam_to_fastq:
        run: ../subworkflows/bam_to_trimmed_fastq.cwl
        in:
            unaligned:
                source: xenosplit/graft_bam
                valueFrom: |
                    ${
                        return {'sequence': {'bam': self} };
                    }
            adapters: trimming_adapters
            adapter_trim_end: trimming_adapter_trim_end
            adapter_min_overlap: trimming_adapter_min_overlap
            max_uncalled: trimming_max_uncalled
            min_readlength: trimming_min_readlength
        out:
            [fastqs, fastq1, fastq2]
    graftbam_star_align_fusion:
        run: ../tools/star_align_fusion.cwl
        in:
            outsam_attrrg_line: outsam_attrrg_line
            star_genome_dir: graft_star_genome_dir
            outfile_name_prefix: graft_outfile_name_prefix
            gtf_file: graft_gtf_file
            fastq:
                source: graftbam_to_fastq/fastq1
                valueFrom: ${ return [self]; }
            fastq2:
                source: graftbam_to_fastq/fastq2
                valueFrom: ${ return [self]; }
        out:
            [aligned_bam, chim_junc, splice_junction_out,log_final]
    star_fusion_detect:
        run: ../tools/star_fusion_detect.cwl
        in:
            star_fusion_genome_dir: star_fusion_genome_dir
            junction_file: graftbam_star_align_fusion/chim_junc
        out:
            [fusion_predictions,fusion_abridged]
    kallisto:
        run: ../tools/kallisto.cwl
        in:
            kallisto_index: kallisto_index
            strand: strand
            fastqs:
                source: graftbam_to_fastq/fastqs
                valueFrom: ${ return [self]; }
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
            input_bam: xenosplit/graft_bam
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
            bam: index_bam/indexed_bam
            reference_annotation: graft_gtf_file
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
    cgpbigwig_bamcoverage:
        run: ../tools/bam_to_bigwig.cwl
        in:
            bam: index_bam/indexed_bam
            reference: reference
        out:
            [outfile]
