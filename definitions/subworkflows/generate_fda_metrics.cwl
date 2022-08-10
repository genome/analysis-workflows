#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Calculate FDA-requested metrics on all aligned and unaligned sequence files"
requirements:
    - class: MultipleInputFeatureRequirement
    - class: SubworkflowFeatureRequirement
    - class: InlineJavascriptRequirement
    - class: StepInputExpressionRequirement
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

outputs:
    unaligned_normal_dna_fastqc_data:
        type: File[]
        outputSource: unaligned_normal_dna_fastqc/fastqc_all_data
    unaligned_normal_dna_table_metrics:
        type: File
        outputSource: unaligned_normal_dna_metrics/unaligned_stats
    unaligned_normal_dna_md5sums:
        type: File
        outputSource: unaligned_normal_dna_md5/md5sum

    unaligned_tumor_dna_fastqc_data:
        type: File[]
        outputSource: unaligned_tumor_dna_fastqc/fastqc_all_data
    unaligned_tumor_dna_table_metrics:
        type: File
        outputSource: unaligned_tumor_dna_metrics/unaligned_stats
    unaligned_tumor_dna_md5sums:
        type: File
        outputSource: unaligned_tumor_dna_md5/md5sum

    unaligned_tumor_rna_fastqc_data:
        type: File[]
        outputSource: unaligned_tumor_rna_fastqc/fastqc_all_data
    unaligned_tumor_rna_table_metrics:
        type: File
        outputSource: unaligned_tumor_rna_metrics/unaligned_stats
    unaligned_tumor_rna_md5sums:
        type: File
        outputSource: unaligned_tumor_rna_md5/md5sum

    aligned_normal_dna_fastqc_data:
        type: File[]
        outputSource: aligned_normal_dna_fastqc/fastqc_all_data
    aligned_normal_dna_table_metrics:
        type: File
        outputSource: aligned_normal_dna_metrics/aligned_stats
    aligned_normal_dna_md5sums:
        type: File
        outputSource: aligned_normal_dna_md5/md5sum

    aligned_tumor_dna_fastqc_data:
        type: File[]
        outputSource: aligned_tumor_dna_fastqc/fastqc_all_data
    aligned_tumor_dna_table_metrics:
        type: File
        outputSource: aligned_tumor_dna_metrics/aligned_stats
    aligned_tumor_dna_md5sums:
        type: File
        outputSource: aligned_tumor_dna_md5/md5sum

    aligned_tumor_rna_fastqc_data:
        type: File[]
        outputSource: aligned_tumor_rna_fastqc/fastqc_all_data
    aligned_tumor_rna_table_metrics:
        type: File
        outputSource: aligned_tumor_rna_metrics/aligned_stats
    aligned_tumor_rna_md5sums:
        type: File
        outputSource: aligned_tumor_rna_md5/md5sum

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
