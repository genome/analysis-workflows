#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "multiqc"
baseCommand: ["/usr/local/bin/multiqc"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
      ramMin: 4000  
    - class: DockerRequirement
      dockerPull: "sridnona/multiqc:1.9"
    - class: StepInputExpressionRequirement

arguments: [
    "--zip-data-dir",
    "-n", "multiqc_out"
    
]
 
inputs:
  inputfiles_array:
    type: File[]
    inputBinding:
        position: 2
  

outputs:
  multiqc_zip:
    type: File
    outputBinding:
        glob: multiqc_out_data.zip
  multiqc_html:
    type: File
    outputBinding:
        glob: multiqc_out.html
