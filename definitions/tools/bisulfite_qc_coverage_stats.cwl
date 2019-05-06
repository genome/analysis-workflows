#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Biscuit QC: Bisulfite conversion rate."
baseCommand: ["/bin/bash", "/opt/biscuit/scripts/Bisulfite_QC_Coveragestats.sh"]
requirements:
            - class: ResourceRequirement
              coresMin: 1
              ramMin: 16000
              tmpdirMin: 20000
            - class: DockerRequirement
              dockerPull: "mgibio/biscuit:0.3.8.2"

inputs:
      vcf:
        type: File
        inputBinding:
          position: 1

      bam:
        type: File
        inputBinding:
          position: 2

      reference:
        type: string
        inputBinding:
          position: 3

      QCannotation:
        type: File
        inputBinding:
          position: 4

outputs:
       bga_bed:
             type: File
             outputBinding:
              glob: "bga.bed"
       cov_dist:
              type: File
              outputBinding:
                glob: "covdist_table.txt"
       bga_bed_dup:
                type: File
                outputBinding:
                  glob: "bga_dup.bed"
       dup_report:
                type: File
                outputBinding:
                  glob: "dup_report.txt"
       cpg_bed:
             type: File
             outputBinding:
              glob: "cpg.bed"
       cov_dist_cpg:
                 type: File
                 outputBinding:
                  glob: "covdist_cpg_table.txt"
       cpg_dist:
              type: File
              outputBinding:
                glob: "cpg_dist_table.txt"
