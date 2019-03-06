#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Run Cell Ranger mkfastq"

baseCommand: ["/opt/cellranger-3.0.1/cellranger", "mkfastq", "--localmem=64", "--localcores=8", "--id=cellranger_output"]

requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/mgi/cellranger:3.0.1"
    - class: ResourceRequirement
      ramMin: 64000
      coresMin: 8

inputs:
    bcl_directory:
        type: Directory
        inputBinding:
            prefix: --run=
            position: 1
            separate: false
        doc: "Directory of the Illumina BCL run"
    simple_sample_csv:  
        type: File
        inputBinding:
            prefix: --simple-csv=
            position: 3
            separate: false
        doc: "input simple sample CSV used to describe the way to demultiplex the flowcell"

outputs:
    samplesheet_csv:
        type: File
        outputBinding:
            glob: "cellranger_output/outs/input_samplesheet.csv"
    fastq_dir:
        type: Directory
        outputBinding:
            glob: "cellranger_output/outs/fastq_path/"
