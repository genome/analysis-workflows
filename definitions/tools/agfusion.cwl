#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: "v1.0"
label: "A tool that annotates STAR gene fusion predictions"

baseCommand: ["/usr/local/bin/agfusion", "batch"]
requirements:
    - class: ShellCommandRequirement 
    - class: ResourceRequirement
      ramMin: 32000
      coresMin: 4 
    - class: DockerRequirement
      dockerPull: rachelbj/agfusion:1.1

arguments: ["-a", "starfusion", "--middlestar"]

inputs:

  fusion_predictions:
    type: File
    inputBinding:
        prefix: '-f'
        position: 1

  agfusion_database:
    type: File
    default: "/opt/agfusiondb/agfusion.homo_sapiens.87.db"
    inputBinding:
        prefix: '-db'
        position: 2

  coding_effect:
    type: boolean
    default: false
    doc: "Annotate all gene transcripts, not just canonical isoforms"
    inputBinding:
        prefix: '--noncanonical'
        position: 3

  output_dir:
    type: string
    default: "agfusion"
    inputBinding:
        prefix: '-o'
        position: 4


outputs:
    annotated_fusion_predictions:
      type: Directory
      outputBinding:
        glob: $(inputs.output_dir)
