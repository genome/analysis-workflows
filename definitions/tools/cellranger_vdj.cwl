#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Run Cell Ranger V(D)J"

baseCommand: ["/apps/cellranger-6.0.0/cellranger", "vdj"]
arguments: ["--id=$(inputs.sample_name)", "--localcores=$(runtime.cores)", "--localmem=$(runtime.ram/1000)"]

requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/alex.paul/cellranger:6.0.0"
    - class: ResourceRequirement
      ramMin: 64000
      coresMin: 8

inputs:
    fastq_directory:
        type: Directory[]
        inputBinding:
            prefix: --fastqs=
            position: 1
            itemSeparator: ","
            separate: false
        doc: "Array of directories containing fastq files"
    reference:
        type: Directory
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
            glob: "$(inputs.sample_name)/outs/"
