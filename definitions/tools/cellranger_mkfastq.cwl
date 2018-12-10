#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Run Cell Ranger mkfastq"

baseCommand: ["cellranger", "mkfastq", "--localmem=64", "--localcores=8"]

requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/ebelter/cellranger:2.2.0"
    - class: ResourceRequirement
      ramMin: 64000
      coresMin: 8

inputs:
    bcl_directory:
        type: string
        inputBinding:
            prefix: --run=
            position: 1
            separate: false
        doc: "Directory of the Illumina BCL run"
    unique_id:
        type: string
        inputBinding:
            prefix: --id=
            position: 2
            separate: false
        doc: "Unique name for the run"
    simple_sample_csv:  
        type: string
        inputBinding:
            prefix: --simple-csv=
            position: 3
            separate: false
        doc: "input simple sample CSV used to describe the way to demultiplex the flowcell"

outputs:
    samplesheet_csv:
        type: File
        outputBinding:
            glob: "$(inputs.unique_id)/outs/input_samplesheet.csv"
    fastq_dir:
        type: Directory
        outputBinding:
            glob: "$(inputs.unique_id)/outs/fastq_path/"