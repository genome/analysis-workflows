#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Run Cell Ranger ATAC Count"

baseCommand: ["/opt/cellranger-atac-1.0.1/cellranger-atac", "count", "--localmem=64", "--localcores=8"]
arguments: ["--id=$(inputs.sample_name)"]

requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/mgi/cellranger-atac:1.0.1"
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
        type: Directory
        inputBinding:
            prefix: --reference=
            position: 2
            separate: false
        doc: "Transcriptome reference compatible with input species and Cell Ranger"
    sample_name:
        type: string
        inputBinding:
            prefix: --sample=
            position: 3
            separate: false
        doc: "Sample name, must be same as name specified in sample sheet in previous cellranger-atac mkfastq step"

outputs:
    out_dir:
        type: Directory
        outputBinding:
            glob: "$(inputs.sample_name)/outs/"
