#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Calculate FDA-requested metrics on all aligned and unaligned sequence files"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: InlineJavascriptRequirement
    - class: StepInputExpressionRequirement
    - class: ScatterFeatureRequirement
    - class: SchemaDefRequirement
      types:
          - $import: ../types/sequence_data.yml

inputs:
    reference:
        type:
            - string
            - File
    unaligned_normal_dna:
        type: ../types/sequence_data.yml#sequence_data[]
    unaligned_tumor_dna:
        type: ../types/sequence_data.yml#sequence_data[]
    unaligned_tumor_rna:
        type: ../types/sequence_data.yml#sequence_data[]
    aligned_normal_dna:
        type: File
    aligned_tumor_dna:
        type: File
    aligned_tumor_rna:
        type: File

    normal_alignment_summary_metrics:
        type: File
    normal_duplication_metrics:
        type: File
    normal_insert_size_metrics:
        type: File
    normal_hs_metrics:
        type: File
    normal_flagstat:
        type: File
    tumor_alignment_summary_metrics:
        type: File
    tumor_duplication_metrics:
        type: File
    tumor_insert_size_metrics:
        type: File
    tumor_hs_metrics:
        type: File
    tumor_flagstat:
        type: File
    rna_metrics:
        type: File

    reference_genome:
        type: string?
    dna_sequencing_platform:
        type: string?
    dna_sequencing_instrument:
        type: string?
    dna_sequencing_kit:
        type: string?
    dna_sequencing_type:
        type: string?
    dna_single_or_paired_end:
        type: string?
    normal_dna_spike_in_error_rate:
        type: string?
    tumor_dna_spike_in_error_rate:
        type: string?
    normal_dna_total_DNA:
        type: string?
    tumor_dna_total_DNA:
        type: string?
    normal_dna_sample_name:
        type: string?
    tumor_dna_sample_name:
        type: string?
    rna_sequencing_platform:
        type: string?
    rna_sequencing_instrument:
        type: string?
    rna_sequencing_kit:
        type: string?
    rna_sequencing_type:
        type: string?
    rna_single_or_paired_end:
        type: string?
    rna_spike_in_error_rate:
        type: string?
    rna_total_RNA:
        type: string?
    rna_RIN_score:
        type: string?
    rna_freq_normalization_method:
        type: string?
    rna_annotation_file:
        type: string?
    rna_sample_name:
        type: string?


outputs:
    unaligned_normal_dna_fastqc_data:
        type: File[]
        outputSource: unaligned_normal_dna_fastqc/fastqc_all_data
    unaligned_normal_dna_table_metrics:
        type: File[]
        outputSource: unaligned_normal_dna_metrics/unaligned_stats
    unaligned_normal_dna_md5sums:
        type: File
        outputSource: unaligned_normal_dna_md5/md5sum
    unaligned_normal_dna_table1:
        type: File
        outputSource: unaligned_normal_dna_table/table

    unaligned_tumor_dna_fastqc_data:
        type: File[]
        outputSource: unaligned_tumor_dna_fastqc/fastqc_all_data
    unaligned_tumor_dna_table_metrics:
        type: File[]
        outputSource: unaligned_tumor_dna_metrics/unaligned_stats
    unaligned_tumor_dna_md5sums:
        type: File
        outputSource: unaligned_tumor_dna_md5/md5sum
    unaligned_tumor_dna_table1:
        type: File
        outputSource: unaligned_tumor_dna_table/table

    unaligned_tumor_rna_fastqc_data:
        type: File[]
        outputSource: unaligned_tumor_rna_fastqc/fastqc_all_data
    unaligned_tumor_rna_table_metrics:
        type: File[]
        outputSource: unaligned_tumor_rna_metrics/unaligned_stats
    unaligned_tumor_rna_md5sums:
        type: File
        outputSource: unaligned_tumor_rna_md5/md5sum
    unaligned_tumor_rna_table1:
        type: File
        outputSource: unaligned_tumor_rna_table/table

    aligned_normal_dna_fastqc_data:
        type: File[]
        outputSource: aligned_normal_dna_fastqc/fastqc_all_data
    aligned_normal_dna_table_metrics:
        type: File
        outputSource: aligned_normal_dna_metrics/aligned_stats
    aligned_normal_dna_md5sums:
        type: File
        outputSource: aligned_normal_dna_md5/md5sum
    aligned_normal_dna_table2:
        type: File
        outputSource: aligned_normal_dna_table/table

    aligned_tumor_dna_fastqc_data:
        type: File[]
        outputSource: aligned_tumor_dna_fastqc/fastqc_all_data
    aligned_tumor_dna_table_metrics:
        type: File
        outputSource: aligned_tumor_dna_metrics/aligned_stats
    aligned_tumor_dna_md5sums:
        type: File
        outputSource: aligned_tumor_dna_md5/md5sum
    aligned_tumor_dna_table2:
        type: File
        outputSource: aligned_tumor_dna_table/table

    aligned_tumor_rna_fastqc_data:
        type: File[]
        outputSource: aligned_tumor_rna_fastqc/fastqc_all_data
    aligned_tumor_rna_table_metrics:
        type: File
        outputSource: aligned_tumor_rna_metrics/aligned_stats
    aligned_tumor_rna_md5sums:
        type: File
        outputSource: aligned_tumor_rna_md5/md5sum
    aligned_tumor_rna_table3:
        type: File
        outputSource: aligned_tumor_rna_table/table

steps:
    unaligned_normal_dna_fastqc:
        run: ../tools/fastqc.cwl
        in:
            input_files:
                source: unaligned_normal_dna
                valueFrom: |
                    ${
                        var files = [];
                        var i;
                        for (i = 0; i < self.length; i = i + 1) {
                            if (self[i].sequence.hasOwnProperty('bam')) { files.push(self[i].sequence.bam); }
                            if (self[i].sequence.hasOwnProperty('fastq1')) { files.push(self[i].sequence.fastq1); }
                            if (self[i].sequence.hasOwnProperty('fastq2')) { files.push(self[i].sequence.fastq2); }
                        }
                        return files;
                    }
        out:
            [fastqc_all_data]
    unaligned_normal_dna_metrics:
        run: ../tools/unaligned_seq_fda_stats.cwl
        scatter: [unaligned_normal_dna]
        scatterMethod: dotproduct
        in:
            input_files:
                source: unaligned_normal_dna
                valueFrom: |
                    ${
                        var files = [];
                        if (self.sequence.hasOwnProperty('bam')) { files.push(self.sequence.bam); }
                        if (self.sequence.hasOwnProperty('fastq1')) { files.push(self.sequence.fastq1); }
                        if (self.sequence.hasOwnProperty('fastq2')) { files.push(self.sequence.fastq2); }
                        return files;
                    }
            output_name:
                default: "normal_dna"
        out:
            [unaligned_stats]
    unaligned_normal_dna_md5:
        run: ../tools/md5sum.cwl
        in:
            input_files:
                source: unaligned_normal_dna
                valueFrom: |
                    ${
                        var files = [];
                        var i;
                        for (i = 0; i < self.length; i = i + 1) {
                            if (self[i].sequence.hasOwnProperty('bam')) { files.push(self[i].sequence.bam); }
                            if (self[i].sequence.hasOwnProperty('fastq1')) { files.push(self[i].sequence.fastq1); }
                            if (self[i].sequence.hasOwnProperty('fastq2')) { files.push(self[i].sequence.fastq2); }
                        }
                        return files;
                    }
            output_name:
                default: "normal_dna_unaligned"
        out:
            [md5sum]
    unaligned_normal_dna_table:
        run: ../tools/generate_fda_tables.cwl
        in:
            table_file_name:
                default: "unaligned_normal_dna_table1.csv"
            table_num:
                default: "table1"
            md5sum_file: unaligned_normal_dna_md5/md5sum
            fastqc_zips: unaligned_normal_dna_fastqc/fastqc_all_data
            unaligned_metrics: unaligned_normal_dna_metrics/unaligned_stats
            sample_name: normal_dna_sample_name
            sequencing_platform: dna_sequencing_platform
            sequencing_instrument: dna_sequencing_instrument
            sequencing_kit: dna_sequencing_kit
            single_or_paired_end: dna_single_or_paired_end
            sequencing_type: dna_sequencing_type
            spike_in_error_rate: normal_dna_spike_in_error_rate
        out:
            [table]


    unaligned_tumor_dna_fastqc:
        run: ../tools/fastqc.cwl
        in:
            input_files:
                source: unaligned_tumor_dna
                valueFrom: |
                    ${
                        var files = [];
                        var i;
                        for (i = 0; i < self.length; i = i + 1) {
                            if (self[i].sequence.hasOwnProperty('bam')) { files.push(self[i].sequence.bam); }
                            if (self[i].sequence.hasOwnProperty('fastq1')) { files.push(self[i].sequence.fastq1); }
                            if (self[i].sequence.hasOwnProperty('fastq2')) { files.push(self[i].sequence.fastq2); }
                        }
                        return files;
                    }
        out:
            [fastqc_all_data]
    unaligned_tumor_dna_metrics:
        run: ../tools/unaligned_seq_fda_stats.cwl
        scatter: [unaligned_tumor_dna]
        scatterMethod: dotproduct
        in:
            input_files:
                source: unaligned_tumor_dna
                valueFrom: |
                    ${
                        var files = [];
                        if (self.sequence.hasOwnProperty('bam')) { files.push(self.sequence.bam); }
                        if (self.sequence.hasOwnProperty('fastq1')) { files.push(self.sequence.fastq1); }
                        if (self.sequence.hasOwnProperty('fastq2')) { files.push(self.sequence.fastq2); }
                        return files;
                    }
            output_name:
                default: "tumor_dna"
        out:
            [unaligned_stats]
    unaligned_tumor_dna_md5:
        run: ../tools/md5sum.cwl
        in:
            input_files:
                source: unaligned_tumor_dna
                valueFrom: |
                    ${
                        var files = [];
                        var i;
                        for (i = 0; i < self.length; i = i + 1) {
                            if (self[i].sequence.hasOwnProperty('bam')) { files.push(self[i].sequence.bam); }
                            if (self[i].sequence.hasOwnProperty('fastq1')) { files.push(self[i].sequence.fastq1); }
                            if (self[i].sequence.hasOwnProperty('fastq2')) { files.push(self[i].sequence.fastq2); }
                        }
                        return files;
                    }
            output_name:
                default: "tumor_dna_unaligned"
        out:
            [md5sum]
    unaligned_tumor_dna_table:
        run: ../tools/generate_fda_tables.cwl
        in:
            table_file_name:
                default: "unaligned_tumor_dna_table1.csv"
            table_num:
                default: "table1"
            md5sum_file: unaligned_tumor_dna_md5/md5sum
            fastqc_zips: unaligned_tumor_dna_fastqc/fastqc_all_data
            unaligned_metrics: unaligned_tumor_dna_metrics/unaligned_stats
            sample_name: tumor_dna_sample_name
            sequencing_platform: dna_sequencing_platform
            sequencing_instrument: dna_sequencing_instrument
            sequencing_kit: dna_sequencing_kit
            single_or_paired_end: dna_single_or_paired_end
            sequencing_type: dna_sequencing_type
            spike_in_error_rate: tumor_dna_spike_in_error_rate
        out:
            [table]


    unaligned_tumor_rna_fastqc:
        run: ../tools/fastqc.cwl
        in: 
            input_files:
                source: unaligned_tumor_rna
                valueFrom: |
                    ${
                        var files = [];
                        var i;
                        for (i = 0; i < self.length; i = i + 1) {
                            if (self[i].sequence.hasOwnProperty('bam')) { files.push(self[i].sequence.bam); }
                            if (self[i].sequence.hasOwnProperty('fastq1')) { files.push(self[i].sequence.fastq1); }
                            if (self[i].sequence.hasOwnProperty('fastq2')) { files.push(self[i].sequence.fastq2); }
                        }
                        return files;
                    }
        out:
            [fastqc_all_data]
    unaligned_tumor_rna_metrics:
        run: ../tools/unaligned_seq_fda_stats.cwl
        scatter: [unaligned_tumor_rna]
        scatterMethod: dotproduct
        in: 
            input_files:
                source: unaligned_tumor_rna
                valueFrom: |
                    ${
                        var files = [];
                        if (self.sequence.hasOwnProperty('bam')) { files.push(self.sequence.bam); }
                        if (self.sequence.hasOwnProperty('fastq1')) { files.push(self.sequence.fastq1); }
                        if (self.sequence.hasOwnProperty('fastq2')) { files.push(self.sequence.fastq2); }
                        return files;
                    }
            output_name:
                default: "tumor_rna"
        out:
            [unaligned_stats]
    unaligned_tumor_rna_md5:
        run: ../tools/md5sum.cwl
        in: 
            input_files:
                source: unaligned_tumor_rna
                valueFrom: |
                    ${
                        var files = [];
                        var i;
                        for (i = 0; i < self.length; i = i + 1) {
                            if (self[i].sequence.hasOwnProperty('bam')) { files.push(self[i].sequence.bam); }
                            if (self[i].sequence.hasOwnProperty('fastq1')) { files.push(self[i].sequence.fastq1); }
                            if (self[i].sequence.hasOwnProperty('fastq2')) { files.push(self[i].sequence.fastq2); }
                        }
                        return files;
                    }
            output_name:
                default: "tumor_rna_unaligned"
        out:
            [md5sum]
    unaligned_tumor_rna_table:
        run: ../tools/generate_fda_tables.cwl
        in:
            table_file_name:
                default: "unaligned_tumor_rna_table1.csv"
            table_num:
                default: "table1"
            md5sum_file: unaligned_tumor_rna_md5/md5sum
            fastqc_zips: unaligned_tumor_rna_fastqc/fastqc_all_data
            unaligned_metrics: unaligned_tumor_rna_metrics/unaligned_stats
            sample_name: rna_sample_name
            sequencing_platform: rna_sequencing_platform
            sequencing_instrument: rna_sequencing_instrument
            sequencing_kit: rna_sequencing_kit
            single_or_paired_end: rna_single_or_paired_end
            sequencing_type: rna_sequencing_type
            spike_in_error_rate: rna_spike_in_error_rate
        out:
            [table]

    aligned_normal_dna_cram_index:
        run: ../tools/index_cram.cwl
        in:
            cram: aligned_normal_dna
        out:
            [indexed_cram]
    aligned_normal_dna_cram_to_bam:
        run: ../tools/cram_to_bam.cwl
        in:
            reference: reference
            cram: aligned_normal_dna_cram_index/indexed_cram
        out:
            [bam]
    aligned_normal_dna_fastqc:
        run: ../tools/fastqc.cwl
        in:
            input_files:
                source: aligned_normal_dna_cram_to_bam/bam
                valueFrom: ${ return [self]; }
        out:
            [fastqc_all_data]
    aligned_normal_dna_metrics:
        run: ../tools/aligned_seq_fda_stats.cwl
        in:
            reference: reference
            input_files:
                source: aligned_normal_dna_cram_to_bam/bam
                valueFrom: ${ return [self]; }
            output_name:
                default: "normal_dna"
        out:
            [aligned_stats]
    aligned_normal_dna_md5:
        run: ../tools/md5sum.cwl
        in:
            input_files:
                source: aligned_normal_dna_cram_to_bam/bam
                valueFrom: ${ return [self]; }
            output_name:
                default: "normal_dna_aligned"
        out:
            [md5sum]
    aligned_normal_dna_table:
        run: ../tools/generate_fda_tables.cwl
        in:
            table_file_name:
                default: "aligned_normal_dna_table2.csv"
            table_num:
                default: "table2"
            md5sum_file: aligned_normal_dna_md5/md5sum
            fastqc_zips: aligned_normal_dna_fastqc/fastqc_all_data
            aligned_metrics: aligned_normal_dna_metrics/aligned_stats
            alignment_summary_metrics: normal_alignment_summary_metrics
            duplication_metrics: normal_duplication_metrics
            insert_size_metrics: normal_insert_size_metrics
            hs_metrics: normal_hs_metrics
            flagstat: normal_flagstat
            sample_name: normal_dna_sample_name
            sequencing_platform: dna_sequencing_platform
            sequencing_instrument: dna_sequencing_instrument
            sequencing_kit: dna_sequencing_kit
            single_or_paired_end: dna_single_or_paired_end
            source:
                default: "Normal sample"
            total_DNA: normal_dna_total_DNA
            reference_genome: reference_genome
        out:
            [table]


    aligned_tumor_dna_cram_index:
        run: ../tools/index_cram.cwl
        in:
            cram: aligned_tumor_dna
        out:
            [indexed_cram]
    aligned_tumor_dna_cram_to_bam:
        run: ../tools/cram_to_bam.cwl
        in:
            reference: reference
            cram: aligned_tumor_dna_cram_index/indexed_cram
        out:
            [bam]
    aligned_tumor_dna_fastqc:
        run: ../tools/fastqc.cwl
        in:
            input_files:
                source: aligned_tumor_dna_cram_to_bam/bam
                valueFrom: ${ return [self]; }
        out:
            [fastqc_all_data]
    aligned_tumor_dna_metrics:
        run: ../tools/aligned_seq_fda_stats.cwl
        in:
            reference: reference
            input_files:
                source: aligned_tumor_dna_cram_to_bam/bam
                valueFrom: ${ return [self]; }
            output_name:
                default: "tumor_dna"
        out:
            [aligned_stats]
    aligned_tumor_dna_md5:
        run: ../tools/md5sum.cwl
        in:
            input_files:
                source: aligned_tumor_dna_cram_to_bam/bam
                valueFrom: ${ return [self]; }
            output_name:
                default: "tumor_dna_aligned"
        out:
            [md5sum]
    aligned_tumor_dna_table:
        run: ../tools/generate_fda_tables.cwl
        in:
            table_file_name:
                default: "aligned_tumor_dna_table2.csv"
            table_num:
                default: "table2"
            md5sum_file: aligned_tumor_dna_md5/md5sum
            fastqc_zips: aligned_tumor_dna_fastqc/fastqc_all_data
            aligned_metrics: aligned_tumor_dna_metrics/aligned_stats
            alignment_summary_metrics: tumor_alignment_summary_metrics
            duplication_metrics: tumor_duplication_metrics
            insert_size_metrics: tumor_insert_size_metrics
            hs_metrics: tumor_hs_metrics
            flagstat: tumor_flagstat
            sample_name: tumor_dna_sample_name
            sequencing_platform: dna_sequencing_platform
            sequencing_instrument: dna_sequencing_instrument
            sequencing_kit: dna_sequencing_kit
            single_or_paired_end: dna_single_or_paired_end
            source:
                default: "Tumor sample"
            total_DNA: tumor_dna_total_DNA
            reference_genome: reference_genome
        out:
            [table]


    aligned_tumor_rna_fastqc:
        run: ../tools/fastqc.cwl
        in:
            input_files:
                source: aligned_tumor_rna
                valueFrom: ${ return [self]; }
        out:
            [fastqc_all_data]
    aligned_tumor_rna_metrics:
        run: ../tools/aligned_seq_fda_stats.cwl
        in:
            reference: reference
            input_files:
                source: aligned_tumor_rna
                valueFrom: ${ return [self]; }
            output_name:
                default: "tumor_rna"
        out:
            [aligned_stats]
    aligned_tumor_rna_md5:
        run: ../tools/md5sum.cwl
        in:
            input_files:
                source: aligned_tumor_rna
                valueFrom: ${ return [self]; }
            output_name:
                default: "tumor_rna_aligned"
        out:
            [md5sum]
    aligned_tumor_rna_table:
        run: ../tools/generate_fda_tables.cwl
        in:
            table_file_name:
                default: "aligned_tumor_rna_table3.csv"
            table_num:
                default: "table3"
            md5sum_file: aligned_tumor_rna_md5/md5sum
            fastqc_zips: aligned_tumor_rna_fastqc/fastqc_all_data
            aligned_metrics: aligned_tumor_rna_metrics/aligned_stats
            rna_metrics: rna_metrics
            sample_name: rna_sample_name
            sequencing_platform: rna_sequencing_platform
            sequencing_instrument: rna_sequencing_instrument
            sequencing_kit: rna_sequencing_kit
            single_or_paired_end: rna_single_or_paired_end
            source:
                default: "Tumor sample"
            total_RNA: rna_total_RNA
            RIN_score: rna_RIN_score
            freq_normalization_method: rna_freq_normalization_method
            annotation_file: rna_annotation_file
            reference_genome: reference_genome
        out:
            [table]
