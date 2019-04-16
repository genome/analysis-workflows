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
              dockerPull: "mgibio/biscuit:v1.4.1"

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
       BgaBed:
             type: File
             outputBinding:
              glob: "bga.bed"
       covdist:
              type: File
              outputBinding:
                glob: "covdist_table.txt"
       BgaBeddup:
                type: File
                outputBinding:
                  glob: "bga_dup.bed"
       Dupreport:
                type: File
                outputBinding:
                  glob: "dup_report.txt"
       CpGbed:
             type: File
             outputBinding:
              glob: "cpg.bed"
       covdistcpg:
                 type: File
                 outputBinding:
                  glob: "covdist_cpg_table.txt"
       cpgdist:
              type: File
              outputBinding:
                glob: "cpg_dist_table.txt"
