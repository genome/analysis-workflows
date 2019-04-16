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
       baseconvertion:
                   type: File
                   outputBinding:
                    glob: "totalBaseConversionRate.txt"
       readconvertion:
                   type: File
                   outputBinding:
                    glob: "totalReadConversionRate.txt"
       CpHretention:
                   type: File
                   outputBinding:
                    glob: "CpHRetentionByReadPos.txt"
       CpGretention:
                   type: File
                   outputBinding:
                    glob: "CpGRetentionByReadPos.txt"
