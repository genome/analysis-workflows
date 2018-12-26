#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Run Cell Ranger V(D)J"

baseCommand: ["cellranger", "vdj", "--localmem=64", "--localcores=8", "--id=cellranger_output"]

requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/ebelter/cellranger:2.2.0"
    - class: ResourceRequirement
      ramMin: 64000
      coresMin: 8

inputs:
    fastq_directory:
        type: Directory
        inputBinding:
            prefix: --fastqs=
            position: 1
            separate: false
        doc: "Directory containing fastq files"
    reference:
        type: string
        inputBinding:
            prefix: --reference=
            position: 2
            separate: false
        doc: "Transcriptome reference compatible with input species and Cell Ranger VDJ"
    sample_name:
        type: string
        inputBinding:
            prefix: --sample=
            position: 3
            separate: false
        doc: "Sample name, must be same as name specified in sample sheet in previous mkfastq step" 
outputs:
    out_dir:
        type: Directory
        outputBinding:
            glob: "cellranger_output/outs/"
