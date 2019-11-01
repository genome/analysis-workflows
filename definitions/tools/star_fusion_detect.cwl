#!/usr/bin/env cwl-runner


class: CommandLineTool
cwlVersion: "v1.0"
label: "STAR-Fusion identify candidate fusion transcript"

baseCommand: ["/usr/local/src/STAR-Fusion/STAR-Fusion"]
requirements:
    - class: ShellCommandRequirement 
    - class: ResourceRequirement
      ramMin: 64000
      coresMin: 10 
    - class: DockerRequirement
      dockerPull: trinityctat/starfusion:1.8.0

arguments: [
    "--CPU", $(runtime.cores)
]


inputs:
  star_fusion_genome_dir:
    type: Directory
    inputBinding:
        position: 2
        prefix: '--genome_lib_dir'
    doc: '
    fusion genome dir
    '

  junction_file:
    type: File
    inputBinding:
        prefix: '-J'
        position: 3

  fusion_output_dir:
    type: string
    default: "STAR-Fusion_outdir"
    inputBinding:
        prefix: '--output_dir'
        position: 4

  star_path:
    type: string
    default: "/usr/local/bin/STAR"
    inputBinding:
        prefix: '--STAR_PATH'
        position: 5

  fusion_inspect:
    type: string
    default: "validate"
    inputBinding:
        prefix: '--FusionInspector'
        position: 6

outputs:
    fusion_predictions:
      type: File
      outputBinding:
        glob: $(inputs.fusion_output_dir +"/star-fusion.fusion_predictions.tsv")
    fusion_abridged:
      type: File
      outputBinding:
        glob: $(inputs.fusion_output_dir +"/star-fusion.fusion_predictions.abridged.tsv")
