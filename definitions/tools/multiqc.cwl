#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "multiqc"
baseCommand: ["/usr/local/bin/multiqc"]
requirements:
    - class: ResourceRequirement
      ramMin: 4000  
    - class: DockerRequirement
      dockerPull: "mgibio/multiqc:1.9"

arguments: [
    "--zip-data-dir",
    "-n", "multiqc_out"
    
]
 
inputs:
  qc_source_files:
    type: File[]
    doc: | "path to file or dir that is supported here https://multiqc.info/#supported-tools"
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
