#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Biscuit QC: Mapping summary."
baseCommand: ["/bin/bash", "/opt/biscuit/scripts/Bisulfite_QC_mappingsummary.sh"]
requirements:
            - class: ResourceRequirement
              coresMin: 1
              ramMin: 16000
              tmpdirMin: 10000
            - class: DockerRequirement
              dockerPull: "mgibio/biscuit"

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
       Strandtable:
             type: File
             outputBinding:
              glob: "strand_table.txt"
       Mappingquality:
             type: File
             outputBinding:
              glob: "mapq_table.txt"
