#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: "v1.0"
label: "A tool that annotates STAR gene fusion predictions"
requirements:
    - class: ResourceRequirement
      ramMin: 32000
      coresMin: 4 
    - class: DockerRequirement
      dockerPull: agfusion:1.3.1-ensembl-95

baseCommand: ["/usr/local/bin/agfusion", "batch"]
arguments: ["-a", "starfusion", "--middlestar"]
inputs:
  fusion_predictions:
    type: File
    inputBinding:
        prefix: '-f'
        position: 1
  agfusion_database:
    type: File
    inputBinding:
        prefix: '-db'
        position: 2
  annotate_noncanonical:
    type: boolean?
    doc: "Annotate all gene transcripts, not just canonical isoforms"
    inputBinding:
        prefix: '--noncanonical'
        position: 3
  output_dir:
    type: string?
    default: "agfusion_results"
    inputBinding:
        prefix: '-o'
        position: 4
outputs:
    annotated_fusion_predictions:
      type: Directory
      outputBinding:
        glob: $(inputs.output_dir)
