#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: "v1.0"
label: "Report metrics on bam or fastq files using fastqc"
requirements:
    - class: ResourceRequirement
      ramMin: 32000
      coresMin: 2
    - class: DockerRequirement
      dockerPull: mgibio/fastqc:0.11.9

baseCommand: ["/usr/bin/perl", "/usr/local/bin/FastQC/fastqc"]
arguments: ["-q", "-o", "$(runtime.outdir)"]

inputs:
  input_files:
    type: File[]
    inputBinding:
        position: 1

outputs:
    fastqc_all_data:
      type: File[]
      outputBinding:
        glob: "*.zip"
