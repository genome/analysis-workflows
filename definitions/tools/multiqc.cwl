#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "multiqc"
baseCommand: ["/usr/local/bin/multiqc"]
requirements:
    - class: ShellCommandRequirement
    - class: ResourceRequirement
    - class: InlineJavascriptRequirement
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
        seperate: false
        position: 2
        valueFrom: |
          ${
          var r = [];
          for(var i=0; i < self.length; i++) {
            r.push(self[i]);
            }
          return r;
          }
  

outputs:
  multiqc_zip:
    type: File
    outputBinding:
        glob: multiqc_out_data.zip
  multiqc_html:
    type: File
    outputBinding:
        glob: multiqc_out.html
