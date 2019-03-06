#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
label: "Run Cell Ranger Count"

baseCommand: ["/opt/cellranger-3.0.1/cellranger", "count", "--localmem=64", "--localcores=8", "--id=cellranger_output"]

requirements:
    - class: DockerRequirement
      dockerPull: "registry.gsc.wustl.edu/mgi/cellranger:3.0.1"
    - class: ResourceRequirement
      ramMin: 64000
      coresMin: 8

inputs:
    chemistry:
        type: string?
        inputBinding:
            prefix: --chemistry=
            position: 1
            separate: false
        default: "auto"
        doc: "Assay configuration used, default 'auto' should usually work without issue"
    fastq_directory:
        type: Directory
        inputBinding:
           prefix: --fastqs=
           position: 2
           separate: false
        doc: "Directory containing fastq files"
    reference:
        type: Directory
        inputBinding:
            prefix: --transcriptome=
            position: 3
            separate: false
        doc: "Transcriptome reference compatible with input species and Cell Ranger"
    sample_name:
        type: string
        inputBinding:
            prefix: --sample=
            position: 4
            separate: false
        doc: "Sample name, must be same as name specified in sample sheet in previous mkfastq step"
outputs:
    out_dir:
        type: Directory
        outputBinding:
            glob: "cellranger_output/outs/"
