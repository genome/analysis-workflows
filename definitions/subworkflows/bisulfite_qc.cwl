#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
label: "Bisulfite QC tools"
requirements:
    - class: SubworkflowFeatureRequirement

inputs:
      vcf:
          type: File
      bam:
          type: File
      reference:
          type: string
      QCannotation:
           type: File
      output_dir:
           type: string
           default: bisulfite_QC

outputs:
    QC_directory:
        type: Directory
        outputSource: gather_to_sub_directory/gathered_directory

steps:
      bisulfite_conversion:
          run: ../tools/bisulfite_qc_conversion.cwl
          in:
            vcf: vcf
            bam: bam
            reference: reference
            QCannotation: QCannotation
          out:
              [base_convertion, read_convertion, cph_retention, cpg_retention]

      mapping_summary:
          run: ../tools/bisulfite_qc_mapping_summary.cwl
          in:
            vcf: vcf
            bam: bam
            reference: reference
            QCannotation: QCannotation
          out:
             [strand_table, mapping_quality]

      cpg_retention_distribution:
          run: ../tools/bisulfite_qc_cpg_retention_distribution.cwl
          in:
            vcf: vcf
            bam: bam
            reference: reference
            QCannotation: QCannotation
          out:
             [cpg_retention_dist]

      coverage_stats:
          run: ../tools/bisulfite_qc_coverage_stats.cwl
          in:
            vcf: vcf
            bam: bam
            reference: reference
            QCannotation: QCannotation
          out:
             [bga_bed, cov_dist, bga_bed_dup, dup_report, cpg_bed, cov_dist_cpg, cpg_dist]

      gather_to_sub_directory:
          run: ../tools/gather_to_sub_directory.cwl
          in:
            outdir: output_dir
            files: [bisulfite_conversion/base_convertion, bisulfite_conversion/read_convertion, bisulfite_conversion/cph_retention, bisulfite_conversion/cpg_retention, mapping_summary/strand_table, mapping_summary/mapping_quality, cpg_retention_distribution/cpg_retention_dist, coverage_stats/bga_bed, coverage_stats/cov_dist, coverage_stats/bga_bed_dup, coverage_stats/dup_report, coverage_stats/cpg_bed, coverage_stats/cov_dist_cpg, coverage_stats/cpg_dist]
          out:
             [gathered_directory]
