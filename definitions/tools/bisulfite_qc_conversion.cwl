#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Biscuit QC: Bisulfite conversion rate."
baseCommand: ["/bin/bash", "/opt/biscuit/scripts/Bisulfite_QC_bisulfiteconversion.sh"]
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
       base_convertion:
                   type: File
                   outputBinding:
                    glob: "totalBaseConversionRate.txt"
       read_convertion:
                   type: File
                   outputBinding:
                    glob: "totalReadConversionRate.txt"
       cph_retention:
                   type: File
                   outputBinding:
                    glob: "CpHRetentionByReadPos.txt"
       cpg_retention:
                   type: File
                   outputBinding:
                    glob: "CpGRetentionByReadPos.txt"
