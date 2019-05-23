#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Running a script to identify lineage of cells"

baseCommand: ["/usr/local/bin/Rscript", "/opt/CellMatch_Haemopedia.r"]

requirements:
    - class: DockerRequirement
      dockerPull: "mgibio/scrna_lineage_inference:0.2"
    - class: ResourceRequirement
      ramMin: 16000
      coresMin: 1

arguments: ["$(inputs.cellranger_out_dir)/filtered_feature_bc_matrix", "cellmatch"]

inputs:
    sample_name:
        type: string
        inputBinding:
            position: -1
    cellranger_out_dir:
        type: Directory
    lineage_reference_data:
        type: string
        inputBinding:
            position: 1
    lineage_min_cells:
        type: int?
        default: 3
        inputBinding:
            position: 2 
    lineage_min_features:
        type: int?
        default: 10
        inputBinding:
            position: 3

outputs:
    cellmatch_out_dir:
        type: Directory
        outputBinding:
            glob: "cellmatch"
